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
import os
import PyAAVF.parser as parser

# pylint: disable=no-self-use,too-few-public-methods

def fhandle(fname, mode='rt'):
    """Return an open file handle."""
    return open(os.path.join(os.path.dirname(__file__), fname), mode)

class TestAAVFSpecs(object):
    """Test whether the AAVF file can be walked through."""
    def test_aavf_1_0(self):
        """Test with AAVF Version 1.0"""
        reader = parser.Reader(fhandle('sample.aavf'))
        assert reader.metadata['fileformat'] == 'AAVFv1.0'

        # test we can walk the file at least
        for line in reader:

            if line.POS == 103:
                assert not line.is_filtered
            else:
                assert line.is_filtered

class TestInfoOrder(object):
    """Test whether items referenced in INFO metadata are ordered correctly"""
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
        reader = parser.Reader(fhandle('sample.aavf', 'r'))
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

class TestInfoTypeCharacter(object):
    """Perfom tests to make sure INFO section is parser and written correctly"""
    def test_parse(self):
        """Test whether the INFO section can be parsed correctly."""
        reader = parser.Reader(fhandle('sample.aavf'))
        record = next(reader)
        assert record.INFO['RC'] == 'tca'
        # the below two RESERVED_INFO constants in the INFO fields have a
        # number of possible values that varies ior is unbounded. Thus, a list
        # is returned.
        assert record.INFO['AC'] == ['tAa']
        assert record.INFO['ACF'] == [0.0031]

    def test_write(self):
        """Test whether the INFO section can be written correctly."""
        reader = parser.Reader(fhandle('sample.aavf'))
        out = StringIO()
        writer = parser.Writer(out, reader)

        records = list(reader)

        for record in records:
            writer.write_record(record)
        out.seek(0)
        reader2 = parser.Reader(out)

        for left, right in zip(records, reader2):
            assert left.INFO == right.INFO
