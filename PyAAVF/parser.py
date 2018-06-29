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

import collections
import itertools
import re

from PyAAVF.model import _Record
from PyAAVF.model import AAVF


MISSING_VALUE = '.'

# Metadata parsers/constants
RESERVED_INFO = {
    'RC': 'String', 'AC': 'String', 'ACC': 'Integer', 'ACF': 'Float',
}

# Metadata lines which are singular
SINGULAR_METADATA = ['fileformat', 'fileDate', 'reference', 'source']

# Conversion between value in file and Python value
FIELD_COUNTS = {
    MISSING_VALUE: None,  # Unknown number of values
}


_Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source',
                                        'version'])
_Filter = collections.namedtuple('Filter', ['id', 'desc'])


class _aavfMetadataParser(object):
    '''Parse the metadata in the header of a AAVF file.'''
    def __init__(self):
        super(_aavfMetadataParser, self).__init__()
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

    # pylint: disable=no-self-use
    def aavf_field_count(self, num_str):
        """Cast aavf header numbers to integer or None"""
        # pylint: disable=no-else-return
        if num_str is None:
            return None
        elif num_str not in FIELD_COUNTS:
            # Fixed, specified number
            return int(num_str)
        else:
            return FIELD_COUNTS[num_str]

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

    def read_meta(self, meta_string):
        '''Read a meta-information META line.'''
        match = self.meta_pattern.match(meta_string)
        if not match:
            # Spec only allows key=value, but we try to be liberal and
            # interpret anything else as key=none (and all values are parsed
            # as strings).
            return meta_string.lstrip('#'), 'none'
        return (match.group('key'), match.group('val'))


# pylint: disable=too-many-instance-attributes,too-many-arguments,too-few-public-methods
class Reader(object):
    """ Reader that can be used for parsing records from AAVF files.
        You must specify file handle."""

    def __init__(self, filehandle):
        """ Create a new Reader for a AAVF file. You must specify file handle.
        """
        super(Reader, self).__init__()

        if not filehandle:
            raise Exception('You must provide a file handle.')

        if filehandle:
            self._reader = filehandle

        self._separator = '\t| +'

        self._row_pattern = re.compile(self._separator)

        self.reader = (line.strip() for line in self._reader if line.strip())

        #: metadata fields from header (string or hash, depending)
        self.metadata = {}
        #: INFO fields from header
        self.infos = {}
        #: FILTER fields from header
        self.filters = {}
        self._header_lines = []
        self.column_headers = []
        self.header_lines = []

    # pylint: disable=too-many-locals
    def read_records(self):
        """ Parse records from a AAVF file, returns an iterable AAVF object which can
            be used to iterate over AAVF records read from a file. The AAVF object
            returned also contains the metadata parsed from the file."""

        self._parse_metadata()
        record_list = []
        line = next(self.reader)

        # add records to the AAVF object until we have reached StopIteration
        while line:
            row = self._row_pattern.split(line.rstrip())
            chrom = row[0]

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
            record_list.append(record)
            try:
                line = next(self.reader)
            except StopIteration:
                line = None

        aavf = AAVF(self.metadata, self.infos, self.filters, self.column_headers, record_list)

        return aavf

    # pylint: disable=dangerous-default-value,no-self-use
    def _map(self, func, iterable, bad=[MISSING_VALUE, '']):
        '''``map``, but make bad values None.'''
        return [func(x) if x not in bad else None
                for x in iterable]

    # pylint: disable=no-self-use,no-else-return
    def _parse_filter(self, filt_str):
        '''Parse the FILTER field of a AAVF entry into a Python list
        '''
        if filt_str == MISSING_VALUE:  # if set to the missing value
            return None
        elif filt_str == 'PASS':
            return []
        else:
            return filt_str.split(';')

    # pylint: disable=too-many-branches
    def _parse_info(self, info_str):
        '''Parse the INFO field of a AAVF entry into a dictionary of Python
        types.
        '''
        if info_str == MISSING_VALUE:  # if set to the missing value
            return {}

        entries = info_str.split(';')
        retdict = {}

        for entry in entries:
            entry = entry.split('=', 1)
            info_id = entry[0]
            try:
                entry_type = self.infos[info_id].type
            except KeyError:
                try:
                    entry_type = RESERVED_INFO[info_id]
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
                if self.infos[info_id].num == 1 and entry_type not in ('Flag', ):
                    val = val[0]
            except KeyError:
                pass

            retdict[info_id] = val

        return retdict

    def _parse_metadata(self):
        '''Parse the information stored in the metainfo of the AAVF and return
           a new iterable AAVF object with the metadata stored in that object.
           Throws an exception if unable to parse header lines of the file.'''

        self.metadata = {}
        self.infos = {}
        self.filters = {}
        self.header_lines = []

        parser = _aavfMetadataParser()
        line = next(self.reader)

        while line.startswith('##'):
            self.header_lines.append(line)

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

        if not line:
            raise Exception("Unable to parse header line in AAVF file.")
        else:
            # pylint: disable=dangerous-default-value
            fields = self._row_pattern.split(line[1:])
            self.column_headers = []
            for field in fields[:9]:
                self.column_headers.append(field.strip())


class Writer(object):
    """Writer for AAVF file. You must supply an output stream such as StringIO,
       and an Reader object to use as a template for the AAVF metadata and
       header. Optionally specify the line terminator."""

    # Reverse keys and values in header field count dictionary
    counts = dict((v, k) for k, v in FIELD_COUNTS.items())

    def __init__(self, stream, template):
        self.template = template
        self.stream = stream

        # Order keys for INFO fields defined in the header (undefined fields
        # get a maximum key).
        self.info_order = collections.defaultdict(
            lambda: len(template.infos),
            dict(zip(template.infos.keys(), itertools.count())))

        two = '##{key}=<ID={0},Description="{1}">\n'
        four = '##{key}=<ID={0},Number={num},Type={2},Description="{3}">\n'
        _num = self._fix_field_count
        for (key, vals) in template.metadata.items():
            if key in SINGULAR_METADATA:
                vals = [vals]
            for val in vals:
                if isinstance(val, dict):
                    values = ','.join('{0}={1}'.format(key, value)
                                      for key, value in val.items())
                    stream.write('##{0}=<{1}>\n'.format(key, values))
                else:
                    stream.write('##{0}={1}\n'.format(key, val))
        for line in template.infos.values():
            stream.write(four.format(key="INFO", *line, num=_num(line.num)))
        for line in template.filters.values():
            stream.write(two.format(key="FILTER", *line))

        self._write_header()

    def _write_header(self):
        """Write the header line of the AAVF output"""
        self.stream.write('#' + '\t'.join(self.template.column_headers) + '\n')

    def write_record(self, record):
        """Write the record into the next line of the AAVF output
           write a record to the file """
        ffs = [record.CHROM, record.GENE, str(record.POS), record.REF]
        ffs += [self._format_alt(record.ALT), self._format_filter(record.FILTER)]
        ffs += [str(record.ALT_FREQ), str(record.COVERAGE), self._format_info(record.INFO)]
        ffs = map(str, ffs)
        self.stream.write('\t'.join(ffs)+'\n')

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

    # pylint: disable=no-else-return
    def _fix_field_count(self, num_str):
        """Restore header number to original state"""
        if num_str not in self.counts:
            return num_str
        else:
            return self.counts[num_str]

    def _format_alt(self, alt):
        """Format ALT data field"""
        return ','.join(self._map(str, alt))

    def _format_filter(self, flt):
        """Get the correctly formatted FILTER data field"""
        if flt == []:
            return 'PASS'
        return self._stringify(flt, none='.', delim=';')

    def _format_info(self, info):
        """Get the correctly formatted INFO data field"""
        if not info:
            return MISSING_VALUE
        return ';'.join(self._stringify_pair(f, info[f]) for f in
                        sorted(info, key=self.order_key))

    def order_key(self, field):
        '''Order by header definition first, alphabetically second.'''
        return self.info_order[field], field

    def _stringify(self, x_var, none='.', delim=','):
        """Convert an object to a string, accounting for missing values"""
        if isinstance(x_var, list):
            return delim.join(self._map(str, x_var, none))
        return str(x_var) if x_var is not None else none

    def _stringify_pair(self, x_var, y_var, none='.', delim=','):
        """Convert a pair of objects to a string (e.g. "X : Y"), accounting
           for missing values."""
        if isinstance(y_var, bool):
            return str(x_var) if y_var else ""
        return "%s=%s" % (str(x_var),
                          self._stringify(y_var, none=none, delim=delim))

    # pylint: disable=no-self-use
    def _map(self, func, iterable, none=MISSING_VALUE):
        '''``map``, but make None values none.'''
        return [func(x_var) if x_var is not None else none
                for x_var in iterable]


def __update_readme():
    import PyAAVF
    open('README.rst', 'w').write(PyAAVF.__doc__)
