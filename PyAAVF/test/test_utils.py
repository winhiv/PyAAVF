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

import os
import PyAAVF.parser as parser
import PyAAVF.utils as utils

# pylint: disable=no-self-use
class TestUtils(object):
    """
    Utilities for AAVF files.
    """

    def fhandle(self, fname, mode='rt'):
        """Return an open file handle."""
        return open(os.path.join(os.path.dirname(__file__), fname), mode)

    def test_walk(self):
        """
        Walk through three readers simultaneously and make sure that the
        output is identical.
        """
        # easy case: all same sites
        reader1 = parser.Reader(self.fhandle('sample.aavf'))
        reader2 = parser.Reader(self.fhandle('sample.aavf'))
        reader3 = parser.Reader(self.fhandle('sample.aavf'))

        number = 0
        for trio in utils.walk_together(reader1, reader2, reader3):
            assert len(trio) == 3
            assert trio[0] == trio[1]
            assert trio[1] == trio[2]
            number += 1
        assert number == 7
