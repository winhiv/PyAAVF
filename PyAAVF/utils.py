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

# Utilities for AAVF files.


def walk_together(*readers):
    """
    Simultaneously iteratate over two or more AAVF readers. For each
    genomic position with a variant, return a list of size equal to the number
    of AAVF readers. This list contains the AAVF record from readers that have
    this variant, and None for readers that don't have it.
    The caller must make sure that inputs are sorted in the same way and use the
    same reference otherwise behaviour is undefined.
    Two variants are considered to be equal if they have the same CHROM, GENE,
    POS, REF, and ALT values.
    """

    nexts = []
    for reader in readers:
        try:
            nexts.append(next(reader))
        except StopIteration:
            nexts.append(None)

    min_k = (None,)   # keep track of the previous min key's contig

    # loop for each record being read by a reader
    while any([reader is not None for reader in nexts]):
        next_index_to_key = dict(
            (index, (r.CHROM, r.GENE, r.POS)) for index, r in enumerate(nexts) if r is not None)

        # key below refers to the values returned by _get_key()
        # only add tuple to list if its key's CHROM values is the same as the
        # CHROM value for the previous min key's contig
        keys_with_prev_contig = [
            key for key in next_index_to_key.values() if key[0] == min_k[0]]

        if any(keys_with_prev_contig):
            min_k = min(keys_with_prev_contig)   # finish previous contig
        else:
            min_k = min(next_index_to_key.values())   # move on to next contig

        # pylint: disable=consider-using-set-comprehension
        min_k_idxs = set([i for i, k in next_index_to_key.items() if k == min_k])
        yield [nexts[i] if i in min_k_idxs else None for i in range(len(nexts))]

        for i in min_k_idxs:
            try:
                nexts[i] = next(readers[i])
            except StopIteration:
                nexts[i] = None
