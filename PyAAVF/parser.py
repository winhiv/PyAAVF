"""
Copyright Government of Canada 2018

Modified by: Matthew Fogel, National Microbiology Laboratory,
             Public Health Agency of Canada

Originally written by: Population Genetics Technologies Ltd,
                       Copyright (c) 2011-2012

Originally written by: John Dougherty,
                       Copyright (c) 2011

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import codecs
import collections
import csv
import gzip
import itertools
import os
import re
import sys

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

try:
    import pysam
except ImportError:
    pysam = None

try:
    import cparse
except ImportError:
    cparse = None

from model import _Record


# Metadata parsers/constants
RESERVED_INFO = {
    'RC': 'String', 'AC': 'String', 'ACC': 'Float', 'ACF': 'Float',
}

# Spec is a bit weak on which metadata lines are singular, like fileformat
# and which can have repeats, like contig
SINGULAR_METADATA = ['fileformat', 'fileDate', 'reference']

# Conversion between value in file and Python value
field_counts = {
    '.': None,  # Unknown number of values
}


_Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source',
                               'version'])
_Filter = collections.namedtuple('Filter', ['id', 'desc'])
_Alt = collections.namedtuple('Alt', ['id', 'desc'])


class _aavf_metadata_parser(object):
    '''Parse the metadata in the header of a AAVF file.'''
    def __init__(self):
        super(_aavf_metadata_parser, self).__init__()
        self.info_pattern = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),\s*
            Number=(?P<number>-?\d+|\.)?,\s*
            Type=(?P<type>Integer|Float|Flag|Character|String),\s*
            Description="(?P<desc>[^"]*)"
            (?:,\s*Source="(?P<source>[^"]*)")?
            (?:,\s*Version="?(?P<version>[^"]*)"?)?
            >''', re.VERBOSE)
        self.filter_pattern = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),\s*
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.meta_pattern = re.compile(r'''##(?P<key>.+?)=(?P<val>.+)''')

    def aavf_field_count(self, num_str):
        """Cast aavf header numbers to integer or None"""
        if num_str is None:
            return None
        elif num_str not in field_counts:
            # Fixed, specified number
            return int(num_str)
        else:
            return field_counts[num_str]

    def read_info(self, info_string):
        '''Read a meta-information INFO line.'''
        match = self.info_pattern.match(info_string)
        if not match:
            raise SyntaxError(
                "One of the INFO lines is malformed: %s" % info_string)

        num = self.aavf_field_count(match.group('number'))

        info = _Info(match.group('id'), num,
                     match.group('type'), match.group('desc'),
                     match.group('source'), match.group('version'))

        return (match.group('id'), info)

    def read_filter(self, filter_string):
        '''Read a meta-information FILTER line.'''
        match = self.filter_pattern.match(filter_string)
        if not match:
            raise SyntaxError(
                "One of the FILTER lines is malformed: %s" % filter_string)

        filt = _Filter(match.group('id'), match.group('desc'))

        return (match.group('id'), filt)

    def read_meta_hash(self, meta_string):
        # assert re.match("##.+=<", meta_string)
        items = meta_string.split('=', 1)
        # Removing initial hash marks
        key = items[0].lstrip('#')
        # N.B., items can have quoted values, so cannot just split on comma
        val = OrderedDict()
        state = 0
        k = ''
        v = ''
        for c in items[1].strip('[<>]'):

            if state == 0:  # reading item key
                if c == '=':
                    state = 1  # end of key, start reading value
                else:
                    k += c  # extend key
            elif state == 1:  # reading item value
                if v == '' and c == '"':
                    v += c  # include quote mark in value
                    state = 2  # start reading quoted value
                elif c == ',':
                    val[k] = v  # store parsed item
                    state = 0  # read next key
                    k = ''
                    v = ''
                else:
                    v += c
            elif state == 2:  # reading quoted item value
                if c == '"':
                    v += c  # include quote mark in value
                    state = 1  # end quoting
                else:
                    v += c
        if k != '':
            val[k] = v
        return key, val

    def read_meta(self, meta_string):
        if re.match("##.+=<", meta_string):
            return self.read_meta_hash(meta_string)
        match = self.meta_pattern.match(meta_string)
        if not match:
            # Spec only allows key=value, but we try to be liberal and
            # interpret anything else as key=none (and all values are parsed
            # as strings).
            return meta_string.lstrip('#'), 'none'
        return match.group('key'), match.group('val')


class Reader(object):
    """ Reader for a AACF file, an iterator returning ``_Record objects`` """

    def __init__(self, fsock=None, filename=None, compressed=None,
                 prepend_chr=False, strict_whitespace=False, encoding='ascii'):
        """ Create a new Reader for a AACF file.
            You must specify either fsock (stream) or filename.  Gzipped
            streams or files are attempted to be recogized by the file
            extension, or gzipped can be forced with ``compressed=True``
            'prepend_chr=True' will put 'chr' before all the CHROM values,
            useful for different sources.
            'strict_whitespace=True' will split records on tabs only (as with
            AACF spec) which allows you to parse files with spaces in the
            sample names.
        """
        super(Reader, self).__init__()

        if not (fsock or filename):
            raise Exception('You must provide at least fsock or filename')

        if fsock:
            self._reader = fsock
            if filename is None and hasattr(fsock, 'name'):
                filename = fsock.name
                if compressed is None:
                    compressed = filename.endswith('.gz')
        elif filename:
            if compressed is None:
                compressed = filename.endswith('.gz')
            self._reader = open(filename, 'rb' if compressed else 'rt')
        self.filename = filename
        if compressed:
            self._reader = gzip.GzipFile(fileobj=self._reader)
            if sys.version > '3':
                self._reader = codecs.getreader(encoding)(self._reader)

        if strict_whitespace:
            self._separator = '\t'
        else:
            self._separator = '\t| +'

        self._row_pattern = re.compile(self._separator)
        self._alt_pattern = re.compile('[\[\]]')

        self.reader = (line.strip() for line in self._reader if line.strip())

        #: metadata fields from header (string or hash, depending)
        self.metadata = None
        #: INFO fields from header
        self.infos = None
        #: FILTER fields from header
        self.filters = None
        self._header_lines = []
        self._column_headers = []
        self._tabix = None
        self._prepend_chr = prepend_chr
        self._parse_metainfo()
        self._format_cache = {}
        self.encoding = encoding

    def __iter__(self):
        return self

    def _parse_metainfo(self):
        '''Parse the information stored in the metainfo of the VCF.
        The end user shouldn't have to use this.  She can access the metainfo
        directly with ``self.metadata``.'''
        for attr in ('metadata', 'infos', 'filters'):
            setattr(self, attr, OrderedDict())

        parser = _aavf_metadata_parser()

        line = next(self.reader)
        while line.startswith('##'):
            self._header_lines.append(line)

            if line.startswith('##INFO'):
                key, val = parser.read_info(line)
                self.infos[key] = val

            elif line.startswith('##FILTER'):
                key, val = parser.read_filter(line)
                self.filters[key] = val

            else:
                key, val = parser.read_meta(line)
                if key in SINGULAR_METADATA:
                    self.metadata[key] = val
                else:
                    if key not in self.metadata:
                        self.metadata[key] = []
                    self.metadata[key].append(val)

            line = next(self.reader)

        fields = self._row_pattern.split(line[1:])
        self._column_headers = fields[:9]

    def _map(self, func, iterable, bad=['.', '']):
        '''``map``, but make bad values None.'''
        return [func(x) if x not in bad else None
                for x in iterable]

    def _parse_filter(self, filt_str):
        '''Parse the FILTER field of a VCF entry into a Python list
        NOTE: this method has a cython equivalent and care must be taken
        to keep the two methods equivalent
        '''
        if filt_str == '.':
            return None
        elif filt_str == 'PASS':
            return []
        else:
            return filt_str.split(';')

    def _parse_info(self, info_str):
        '''Parse the INFO field of a VCF entry into a dictionary of Python
        types.
        '''
        if info_str == '.':
            return {}

        entries = info_str.split(';')
        retdict = {}

        for entry in entries:
            entry = entry.split('=', 1)
            ID = entry[0]
            try:
                entry_type = self.infos[ID].type
            except KeyError:
                try:
                    entry_type = RESERVED_INFO[ID]
                except KeyError:
                    if entry[1:]:
                        entry_type = 'String'
                    else:
                        entry_type = 'Flag'

            if entry_type == 'Integer':
                vals = entry[1].split(',')
                try:
                    val = self._map(int, vals)
                # Allow specified integers to be flexibly parsed as floats.
                # Handles cases with incorrectly specified header types.
                except ValueError:
                    val = self._map(float, vals)
            elif entry_type == 'Float':
                vals = entry[1].split(',')
                val = self._map(float, vals)
            elif entry_type == 'Flag':
                val = True
            elif entry_type in ('String', 'Character'):
                try:
                    vals = entry[1].split(',')  # commas are reserved
                    # characters indicating multiple values
                    val = self._map(str, vals)
                except IndexError:
                    entry_type = 'Flag'
                    val = True

            try:
                if self.infos[ID].num == 1 and entry_type not in ('Flag', ):
                    val = val[0]
            except KeyError:
                pass

            retdict[ID] = val

        return retdict

    def next(self):
        '''Return the next record in the file.'''
        line = next(self.reader)
        row = self._row_pattern.split(line.rstrip())
        chrom = row[0]
        if self._prepend_chr:
            chrom = 'chr' + chrom

        gene = str(row[1])

        pos = int(row[2])

        ref = row[3]
        alt = row[4].split(',')

        filt = self._parse_filter(row[5])

        alt_freq = float(row[6])

        coverage = int(row[7])

        info = self._parse_info(row[8])

        record = _Record(chrom, gene, pos, ref, alt, filt, alt_freq, coverage,
                         info)

        return record

    def fetch(self, chrom, start=None, end=None):
        """ Fetches records from a tabix-indexed AAVF file and returns an
            iterable of ``_Record`` instances
            chrom must be specified.
            The start and end coordinates are in the zero-based,
            half-open coordinate system, similar to ``_Record.start`` and
            ``_Record.end``. The very first base of a chromosome is
            index 0, and the the region includes bases up to, but not
            including the base at the end coordinate. For example
            ``fetch('4', 10, 20)`` would include all variants
            overlapping a 10 base pair region from the 11th base of
            through the 20th base (which is at index 19) of chromosome
            4. It would not include the 21st base (at index 20). See
            http://genomewiki.ucsc.edu/index.php/Coordinate_Transforms
            for more information on the zero-based, half-open coordinate
            system.
            If end is omitted, all variants from start until the end of
            the chromosome chrom will be included.
            If start and end are omitted, all variants on chrom will be
            returned.
            requires pysam
        """
        if not pysam:
            raise Exception('pysam not available, try "pip install pysam"?')
        if not self.filename:
            raise Exception('Please provide a filename (or a "normal" fsock)')

        if not self._tabix:
            self._tabix = pysam.Tabixfile(self.filename,
                                          encoding=self.encoding)

        if self._prepend_chr and chrom[:3] == 'chr':
            chrom = chrom[3:]

        self.reader = self._tabix.fetch(chrom, start, end)
        return self


class Writer(object):
    """AAVF Writer. On Windows Python 2, open stream with 'wb'."""

    # Reverse keys and values in header field count dictionary
    counts = dict((v, k) for k, v in field_counts.iteritems())

    def __init__(self, stream, template, lineterminator="\n"):
        self.writer = csv.writer(stream, delimiter="\t",
                                 lineterminator=lineterminator,
                                 quotechar='', quoting=csv.QUOTE_NONE)
        self.template = template
        self.stream = stream

        # Order keys for INFO fields defined in the header (undefined fields
        # get a maximum key).
        self.info_order = collections.defaultdict(
            lambda: len(template.infos),
            dict(zip(template.infos.iterkeys(), itertools.count())))

        two = '##{key}=<ID={0},Description="{1}">\n'
        four = '##{key}=<ID={0},Number={num},Type={2},Description="{3}">\n'
        _num = self._fix_field_count
        for (key, vals) in template.metadata.iteritems():
            if key in SINGULAR_METADATA:
                vals = [vals]
            for val in vals:
                if isinstance(val, dict):
                    values = ','.join('{0}={1}'.format(key, value)
                                      for key, value in val.items())
                    stream.write('##{0}=<{1}>\n'.format(key, values))
                else:
                    stream.write('##{0}={1}\n'.format(key, val))
        for line in template.infos.itervalues():
            stream.write(four.format(key="INFO", *line, num=_num(line.num)))
        for line in template.filters.itervalues():
            stream.write(two.format(key="FILTER", *line))

        self._write_header()

    def _write_header(self):
        # TODO: write INFO, etc
        self.stream.write('#' + '\t'.join(self.template._column_headers
                                          + self.template.samples) + '\n')

    def write_record(self, record):
        """ write a record to the file """
        ffs = self._map(str, [record.CHROM, record.GENE, record.POS]) \
            + [record.REF, record.ALT, self._format_filter(record.FILTER),
               record.ALT_FREQ, record.COVERAGE,
               self._format_info(record.INFO)]

        self.writer.writerow(ffs)

    def flush(self):
        """Flush the writer"""
        try:
            self.stream.flush()
        except AttributeError:
            pass

    def close(self):
        """Close the writer"""
        try:
            self.stream.close()
        except AttributeError:
            pass

    def _fix_field_count(self, num_str):
        """Restore header number to original state"""
        if num_str not in self.counts:
            return num_str
        else:
            return self.counts[num_str]

    def _format_alt(self, alt):
        return ','.join(self._map(str, alt))

    def _format_filter(self, flt):
        if flt == []:
            return 'PASS'
        return self._stringify(flt, none='.', delim=';')

    def _format_info(self, info):
        if not info:
            return '.'

        def order_key(field):
            # Order by header definition first, alphabetically second.
            return self.info_order[field], field
        return ';'.join(self._stringify_pair(f, info[f]) for f in
                        sorted(info, key=order_key))

    def _format_sample(self, fmt, sample):
        if hasattr(sample.data, 'GT'):
            gt = sample.data.GT
        else:
            gt = './.' if 'GT' in fmt else ''

        result = [gt] if gt else []
        # Following the VCF spec, GT is always the first item whenever it is
        # present.
        for field in sample.data._fields:
            value = getattr(sample.data, field)
            if field == 'GT':
                continue
            if field == 'FT':
                result.append(self._format_filter(value))
            else:
                result.append(self._stringify(value))
        return ':'.join(result)

    def _stringify(self, x, none='.', delim=','):
        if type(x).isinstance(type([])):
            return delim.join(self._map(str, x, none))
        return str(x) if x is not None else none

    def _stringify_pair(self, x, y, none='.', delim=','):
        if isinstance(y, bool):
            return str(x) if y else ""
        return "%s=%s" % (str(x), self._stringify(y, none=none, delim=delim))

    def _map(self, func, iterable, none='.'):
        '''``map``, but make None values none.'''
        return [func(x) if x is not None else none
                for x in iterable]


def __update_readme():
    import vcf
    file('README.rst', 'w').write(vcf.__doc__)
