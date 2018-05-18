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


class _Record(object):
    """ Equivalent to a row in an AAVF file.
        The standard AAVF fields CHROM, GENE, POS, REF, ALT, FILTER, ALT_FREQ,
        COVERAGE and INFO are available as properties.
    """

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

    # For Python 2
    def __cmp__(self, other):
        """__cmp___"""
        return cmp((self.CHROM, self.POS),
                   (getattr(other, "CHROM", None),
                    getattr(other, "POS", None)))

    # For Python 3
    def __eq__(self, other):
        """ _Records are equal if they describe the same variant (same position, amino acids) """
        return (self.CHROM == getattr(other, "CHROM", None) and
                self.POS == getattr(other, "POS", None) and
                self.REF == getattr(other, "REF", None) and
                self.ALT == getattr(other, "ALT", None))

    # For Python 3
    def __lt__(self, other):
        """__lt__"""
        lhs = (self.CHROM, self.POS)
        rhs = (getattr(other, "CHROM", None), getattr(other, "POS", None))
        return lhs < rhs

    def __str__(self):
        """str"""
        return "Record(CHROM=%(CHROM)s, GENE=%(GENE)s, POS=%(POS)s, REF=%(REF)s" % self.__dict__

    def add_filter(self, flt):
        """add filter"""
        if self.FILTER is None:
            self.FILTER = [flt]
        else:
            self.FILTER.append(flt)

    def add_info(self, info, value=True):
        """add_info"""
        self.INFO[info] = value
