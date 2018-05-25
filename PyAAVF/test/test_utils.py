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
"""
Utilities for AAVF files.
"""

from StringIO import StringIO
import PyAAVF.parser as parser
import PyAAVF.utils as utils
import os

class TestUtils():

    def fh(self, fname, mode='rt'):
        return open(os.path.join(os.path.dirname(__file__), fname), mode)

    def test_walk(self):
        # easy case: all same sites
        reader1 = parser.Reader(self.fh('sample.aavf'))
        reader2 = parser.Reader(self.fh('sample.aavf'))
        reader3 = parser.Reader(self.fh('sample.aavf'))

        n = 0
        for x in utils.walk_together(reader1, reader2, reader3):
            assert len(x) == 3
            assert x[0] == x[1]
            assert x[1] == x[2]
            n += 1
        assert n == 7
