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


class AAVF(object):
    '''An iterable AAVF object that contains the metadata from an AAVF file
       and an list of AAVF records.'''

    # pylint: disable=too-many-arguments,too-few-public-methods
    def __init__(self, metadata, infos, filters, column_headers, records):
        """ Create a new AAVF objeect, an iterator of _Record objects.
        """
        #: METADATA fields from header
        self.metadata = metadata
        #: INFO fields from header
        self.infos = infos
        #: FILTER fields from header
        self.filters = filters

        self.column_headers = column_headers

        # a list of _Record objects
        self.records = records

        self.record_len = len(self.records)

        self.index = -1

    def __iter__(self):
        return self

    def __next__(self):
        if self.index == (self.record_len - 1):
            raise StopIteration
        self.index = self.index + 1
        return self.records[self.index]

    def next(self):
        """Get next item in iterable."""
        return self.__next__()


class _Record(object):
    """ Equivalent to a row in an AAVF file.
        The standard AAVF fields CHROM, GENE, POS, REF, ALT, FILTER, ALT_FREQ,
        COVERAGE and INFO are available as properties.
    """
    # pylint: disable=invalid-name,too-many-instance-attributes,too-many-arguments
    def __init__(self, CHROM, GENE, POS, REF, ALT, FILTER, ALT_FREQ, COVERAGE,
                 INFO):
        """init"""
        self.CHROM = CHROM
        self.GENE = GENE
        self.POS = POS
        self.REF = REF
        self.ALT = ALT
        self.FILTER = FILTER
        self.ALT_FREQ = ALT_FREQ
        self.COVERAGE = COVERAGE
        self.INFO = INFO

    # For Python 3
    def __eq__(self, other):
        """ _Records are equal if they describe the same variant (same position, amino acids).
            If AC is present in the INFO column of both _Records then INFO[AC] must be equal
            for _Records to be equal.
        """
        ac_present = False
        if "AC" in self.INFO and "AC" not in other.INFO:
            return False
        if "AC" not in self.INFO and "AC" in other.INFO:
            return False

        if "AC" in self.INFO and "AC" in other.INFO:
            ac_present = True  # AC is present in both INFO fields

        equal_rec = (self.CHROM == getattr(other, "CHROM", None) and
                     self.GENE == getattr(other, "GENE", None) and
                     self.POS == getattr(other, "POS", None) and
                     self.REF == getattr(other, "REF", None) and
                     self.ALT == getattr(other, "ALT", None))

        if ac_present:
            equal_rec = equal_rec and \
                        (self.INFO["AC"] == getattr(other, "INFO", None).get("AC", None))

        return equal_rec

    # For Python 3
    def __lt__(self, other):
        """__lt__"""
        lhs = (self.CHROM, self.GENE, self.POS)
        rhs = (getattr(other, "CHROM", None),
               getattr(other, "GENE", None),
               getattr(other, "POS", None))
        return lhs < rhs

    def __str__(self):
        """str"""
        return (("Record(CHROM=%(CHROM)s, GENE=%(GENE)s, POS=%(POS)s, REF=%(REF)s, " +
                 "ALT=%(ALT)s, FILTER=%(FILTER)s, ALT_FREQ=%(ALT_FREQ)s, " +
                 "COVERAGE=%(COVERAGE)s, INFO=%(INFO)s)") % self.__dict__)

    def add_filter(self, flt):
        """add filter"""
        if self.FILTER is None:
            self.FILTER = [flt]
        else:
            self.FILTER.append(flt)

    def add_info(self, info, value=True):
        """add_info"""
        self.INFO[info] = value

    # pylint: disable=no-else-return
    @property
    def is_filtered(self):
        """ Return True if a variant has been filtered """
        filt = self.FILTER
        if filt is None or not filt:  # FILTER is not set or set to PASS
            return False
        else:
            return True
