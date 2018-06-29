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

from io import StringIO
import os
import PyAAVF.parser as parser
from PyAAVF.model import AAVF
from PyAAVF.model import _Record

# pylint: disable=no-self-use,too-few-public-methods

def fhandle(fname, mode='r'):
    """Return an open file handle."""
    dirname = os.path.dirname(__file__)
    file_d = open(os.path.join(dirname, fname), mode)
    return file_d

class TestAAVFSpecs(object):
    """Test whether the AAVF file can be walked through."""
    def test_aavf_1_0(self):
        """Test with AAVF Version 1.0"""
        reader = parser.Reader()
        aavf_obj = reader.read_records(fhandle('sample.aavf'))

        assert 'fileformat' in aavf_obj.metadata.keys(), "Metadata should contain fileformat," + \
               "metadata dict is %s" % aavf_obj.metadata.items()

        assert aavf_obj.metadata['fileformat'] == 'AAVFv1.0'

        # test we can walk the file at least
        for line in aavf_obj:

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
        reader = parser.Reader()
        aavf_obj = reader.read_records(fhandle('sample.aavf', 'r'))
        out = StringIO()
        writer = parser.Writer(out, aavf_obj)

        for record in aavf_obj:
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
        reader = parser.Reader()
        aavf_obj = reader.read_records(fhandle('sample.aavf'))
        record = next(aavf_obj)
        assert record.INFO['RC'] == 'tca', "record.INFO['RC'] should be 'tca'" + \
               ", record.INFO is %s" %  record.INFO
        # the below two RESERVED_INFO constants in the INFO fields have a
        # number of possible values that varies ior is unbounded. Thus, a list
        # is returned.
        assert record.INFO['AC'] == ['tAa']
        assert record.INFO['ACF'] == [0.0031]

    def test_write(self):
        """Test whether the INFO section can be written correctly."""
        reader = parser.Reader()
        aavf_obj = reader.read_records(fhandle('sample.aavf'))
        out = StringIO()
        writer = parser.Writer(out, aavf_obj)

        records = list(aavf_obj)

        for record in records:
            writer.write_record(record)
        out.seek(0)
        reader2 = parser.Reader()
        aavf_obj2 = reader2.read_records(out)

        for left, right in zip(records, aavf_obj2):
            assert left.INFO == right.INFO, "left.INFO is %s and right.INFO is %s" \
                   % (left.INFO, right.INFO)

class TestWriter(object):
    """Perfom tests to make sure that the Writer is performing as expected"""
    def test_write_to_file(self):
        """Test whether writes to file work as expected."""
        reader = parser.Reader()
        aavf_obj = reader.read_records(fhandle('sample.aavf'))
        out = fhandle('sampleoutput.aavf', "w+")
        writer = parser.Writer(out, aavf_obj)

        records = list(aavf_obj)

        for record in records:
            writer.write_record(record)

        out.close()
        reader1 = parser.Reader().read_records(fhandle('sample.aavf'))

        reader2 = parser.Reader().read_records(fhandle('sampleoutput.aavf'))
        assert len(list(reader1)) == len(list(reader2))
        # all data lines should be read from the sample file

        reader2 = parser.Reader().read_records(fhandle('sampleoutput.aavf'))
        for left, right in zip(reader1, reader2):
            assert left.INFO == right.INFO

class TestReader(object):
    """Perfom tests to make sure that the Reader is performing as expected"""
    def test_read_from_file(self):
        """Test whether reads from file work as expected and if the AAVF record
           object returned is correct."""
        aavf = parser.Reader().read_records(fhandle('sample.aavf'))
        record_list = [record for record in aavf]

        assert isinstance(aavf, AAVF)

        assert aavf.metadata.get("fileformat") == "AAVFv1.0", \
               "fileformat should be AAVFv1.0, metadata is %s" % aavf.metadata
        assert aavf.metadata.get("fileDate") == "20180501", \
               "filedate should be 20180501, metadata is %s" % aavf.metadata
        assert aavf.metadata.get("source") == "myProgramV1.0", \
               "source should be myProgramV1.0, metadata is %s" % aavf.metadata
        assert aavf.metadata.get("reference") == "hxb2.fas", \
               "reference should be hxb2.fas, metadata is %s" % aavf.metadata
        assert aavf.infos
        assert aavf.filters

        assert len(record_list) == 7
        # all data lines should be the same as in the sample file
        for record in record_list:
            assert isinstance(record, _Record)
