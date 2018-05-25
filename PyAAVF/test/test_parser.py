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

from StringIO import StringIO
import PyAAVF.parser as parser
import os


def fh(fname, mode='rt'):
    return open(os.path.join(os.path.dirname(__file__), fname), mode)

class TestInfoOrder():

    def _assert_order(self, definitions, fields):
        """
        Elements common to both lists should be in the same order. Elements
        only in `fields` should be last and in alphabetical order.
        """
        used_definitions = [d for d in definitions if d in fields]
        assert used_definitions == fields[:len(used_definitions)]
        assert fields[len(used_definitions):] == sorted(fields[len(used_definitions):])

    def test_writer(self):
        """
        Order of INFO fields should be compatible with the order of their
        definition in the header and undefined fields should be last and in
        alphabetical order.
        """
        reader = parser.Reader(fh('sample.aavf', 'r'))
        out = StringIO()
        writer = parser.Writer(out, reader, lineterminator='\n')

        for record in reader:
            writer.write_record(record)
        out.seek(0)
        out_str = out.getvalue()

        definitions = []
        for line in out_str.split('\n'):
            if line.startswith('##INFO='):
                definitions.append(line.split('ID=')[1].split(',')[0])
            if not line or line.startswith('#'):
                continue
            fields = [f.split('=')[0] for f in line.split('\t')[7].split(';')]
            self._assert_order(definitions, fields)
