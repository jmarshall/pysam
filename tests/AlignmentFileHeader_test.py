import os
import re
import copy
from collections import OrderedDict as odict
import pytest

import pysam
import pysam.samtools
from TestUtils import make_data_files, BAM_DATADIR


def setUpModule():
    make_data_files(BAM_DATADIR)


class TestHeaderConstruction:
    """testing header construction."""

    header_dict = odict(
        [('SQ', [odict([('LN', 1575), ('SN', 'chr1'), ('AH', 'chr1:5000000-5010000')]),
                 odict([('LN', 1584), ('SN', 'chr2'), ('AH', '*'), ('TP', 'linear')])]),
         ('RG', [odict([('LB', 'SC_1'), ('ID', 'L1'), ('SM', 'NA12891'),
                        ('PU', 'SC_1_10'), ("CN", "name:with:colon")]),
                 odict([('LB', 'SC_2'), ('ID', 'L2'), ('SM', 'NA12891'),
                        ('PU', 'SC_2_12'), ("CN", "name:with:colon")])]),
         ('PG', [odict([('ID', 'P1'), ('VN', '1.0')]),
                 odict([('ID', 'P2'), ('VN', '1.1')])]),
         ('HD', odict([('VN', '1.0')])),
         ('CO', ['this is a comment', 'this is another comment']),
        ])

    header_text = ("@HD\tVN:1.0\n"
                   "@SQ\tSN:chr1\tLN:1575\tAH:chr1:5000000-5010000\n"
                   "@SQ\tSN:chr2\tLN:1584\tAH:*\tTP:linear\n"
                   "@RG\tID:L1\tPU:SC_1_10\tLB:SC_1\tSM:NA12891\tCN:name:with:colon\n"
                   "@RG\tID:L2\tPU:SC_2_12\tLB:SC_2\tSM:NA12891\tCN:name:with:colon\n"
                   "@PG\tID:P1\tVN:1.0\n"
                   "@PG\tID:P2\tVN:1.1\n"
                   "@CO\tthis is a comment\n"
                   "@CO\tthis is another comment\n")

    header_from_references = odict(
        [('SQ', [odict([('LN', 1575), ('SN', 'chr1')]),
                 odict([('LN', 1584), ('SN', 'chr2')])]),
         ('RG', [odict([('LB', 'SC_1'), ('ID', 'L1'), ('SM', 'NA12891'),
                        ('PU', 'SC_1_10'), ("CN", "name:with:colon")]),
                 odict([('LB', 'SC_2'), ('ID', 'L2'), ('SM', 'NA12891'),
                        ('PU', 'SC_2_12'), ("CN", "name:with:colon")])]),
         ('PG', [odict([('ID', 'P1'), ('VN', '1.0')]),
                 odict([('ID', 'P2'), ('VN', '1.1')])]),
         ('HD', odict([('VN', '1.0')])),
         ('CO', ['this is a comment', 'this is another comment']),
        ])

    header_without_text = odict(
        [('SQ', [odict([('LN', 1575), ('SN', 'chr1')]),
                 odict([('LN', 1584), ('SN', 'chr2')])]),
        ])

    def compare_headers(self, test_header, ref_header=None):
        '''compare two headers a and b.'''
        test_header_dict = test_header.as_dict()
        if ref_header is None:
            ref_header = self.header_dict

        for ak, av in test_header_dict.items():
            assert ak in self.header_dict
            assert av == ref_header[ak]
        for ak, av in ref_header.items():
            assert ak in test_header_dict
            assert av == test_header_dict[ak]

    def check_name_mapping(self, test_header):
        for x, y in enumerate(("chr1", "chr2")):
            tid = test_header.get_tid(y)
            ref = test_header.get_reference_name(x)
            assert tid == x
            assert ref == y

        assert test_header.get_tid("chr?") == -1
        with pytest.raises(ValueError): test_header.get_reference_name(2)

    def test_header_constructed_from_dict(self):
        header = pysam.AlignmentHeader.from_dict(self.header_dict)
        self.compare_headers(header)
        self.check_name_mapping(header)

    def test_header_constructed_from_text(self):
        header = pysam.AlignmentHeader.from_text(self.header_text)
        self.compare_headers(header)
        self.check_name_mapping(header)

    def test_header_constructed_from_header(self):
        header = pysam.AlignmentHeader.from_text(self.header_text)
        self.compare_headers(header.copy())
        self.check_name_mapping(header)

    def test_header_constructed_from_references(self):
        text = re.sub("@SQ[^\n]+\n", "", self.header_text)
        assert "@SQ" not in text
        header = pysam.AlignmentHeader.from_references(
            reference_names=["chr1", "chr2"],
            reference_lengths=[1575, 1584],
            text=text)
        self.compare_headers(header, self.header_from_references)
        self.check_name_mapping(header)

    def test_header_constructed_from_references_without_text(self):
        header = pysam.AlignmentHeader.from_references(
            reference_names=["chr1", "chr2"],
            reference_lengths=[1575, 1584])
        self.compare_headers(header, self.header_without_text)
        self.check_name_mapping(header)


class TestHeaderSAM:
    """testing header manipulation"""

    header = {'SQ': [{'LN': 1575, 'SN': 'chr1', 'AH': 'chr1:5000000-5010000'},
                     {'LN': 1584, 'SN': 'chr2', 'AH': '*'}],
              'RG': [{'LB': 'SC_1', 'ID': 'L1', 'SM': 'NA12891',
                      'PU': 'SC_1_10', "CN": "name:with:colon"},
                     {'LB': 'SC_2', 'ID': 'L2', 'SM': 'NA12891',
                      'PU': 'SC_2_12', "CN": "name:with:colon"}],
              'PG': [{'ID': 'P1', 'VN': '1.0'}, {'ID': 'P2', 'VN': '1.1'}],
              'HD': {'VN': '1.0'},
              'CO': ['this is a comment', 'this is another comment'],
              }

    def compare_headers(self, a, b):
        '''compare two headers a and b.'''
        for ak, av in a.items():
            assert ak in b
            assert av == b[ak]

    def setup_method(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex3.sam"),
            "r")

    def test_header_content_is_as_expected(self):
        self.compare_headers(self.header, self.samfile.header.to_dict())
        self.compare_headers(self.samfile.header.to_dict(), self.header)

    def test_text_access_works(self):
        assert self.samfile.text == self.samfile.header.__str__()

    def test_name_mapping(self):
        for x, y in enumerate(("chr1", "chr2")):
            tid = self.samfile.gettid(y)
            ref = self.samfile.getrname(x)
            assert tid == x
            assert ref == y

        assert self.samfile.gettid("chr?") == -1
        with pytest.raises(ValueError): self.samfile.getrname(2)

    def test_dictionary_access_works(self):
        for key in self.header.keys():
            self.compare_headers({key: self.header[key]},
                                 {key: self.samfile.header[key]})

    def test_dictionary_setting_raises_error(self):
        with pytest.raises(TypeError):
            self.samfile.header.__setitem__("CO", ["This is a final comment"])

    def test_dictionary_len_works(self):
        assert len(self.header) == len(self.samfile.header)

    def test_dictionary_keys_works(self):
        # sort for py2.7
        assert sorted(self.header.keys()) == sorted(self.samfile.header.keys())

    def test_dictionary_values_works(self):
        assert len(self.header.values()) == len(self.samfile.header.values())

    def test_dictionary_get_works(self):
        assert self.header.get("HD") == {'VN': '1.0'}
        assert self.header.get("UK", "xyz") == "xyz"
        assert self.header.get("UK") is None

    def test_dictionary_contains_works(self):
        assert "HD" in self.header
        assert "UK" not in self.header

    def teardown_method(self):
        self.samfile.close()


class TestHeaderBAM(TestHeaderSAM):

    def setup_method(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex3.bam"),
            "rb")


class TestHeaderCRAM(TestHeaderSAM):

    def setup_method(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex3.cram"),
            "rc")

    def compare_headers(self, a, b):
        '''compare two headers a and b.'''
        def _strip(dd):
            for x in dd:
                for y in ("M5", "UR"):
                    if y in x:
                        del x[y]
        for ak, av in a.items():
            _strip(av)
            assert ak in b
            _strip(b[ak])

            assert av == b[ak]


class TestHeaderFromRefs:
    '''see issue 144

    reference names need to be converted to string for python 3
    '''

    # def testHeader(self, tmp_path):
    #     refs = ['chr1', 'chr2']
    #     tmpfile = str(tmp_path / f"tmp_{id(self)}")
    #     s = pysam.AlignmentFile(tmpfile, 'wb',
    #                       referencenames=refs,
    #                       referencelengths=[100]*len(refs))
    #     s.close()
    #
    #     assert checkBinaryEqual('issue144.bam', tmpfile), 'bam files differ'


class TestHeaderWriteRead:
    header = {'SQ': [{'LN': 1575, 'SN': 'chr1'},
                     {'LN': 1584, 'SN': 'chr2'}],
              'RG': [{'LB': 'SC_1', 'ID': 'L1', 'SM': 'NA12891',
                      'PU': 'SC_1_10', "CN": "name:with:colon"},
                     {'LB': 'SC_2', 'ID': 'L2', 'SM': 'NA12891',
                      'PU': 'SC_2_12', "CN": "name:with:colon"}],
              'PG': [{'ID': 'P1', 'VN': '1.0', 'CL': 'tool'},
                     {'ID': 'P2', 'VN': '1.1', 'CL': 'tool --option argument', 'PP': 'P1'}],
              'HD': {'VN': '1.0'},
              'CO': ['this is a comment', 'this is another comment'],
              }

    def compare_headers(self, a, header_b):
        '''compare two headers a and b.

        Ignore M5 and UR field as they are set application specific.
        '''
        b = header_b.to_dict()
        for ak, av in a.items():
            assert ak in b
            assert len(av) == len(b[ak])

            for row_a, row_b in zip(av, b[ak]):
                if isinstance(row_b, dict):
                    for x in ["M5", "UR"]:
                        try:
                            del row_b[x]
                        except KeyError:
                            pass
                assert row_a == row_b

    def check_read_write(self, fn, flag_write, header):
        with pysam.AlignmentFile(
                fn,
                flag_write,
                header=header,
                reference_filename=os.path.join(BAM_DATADIR, "ex1.fa")) as outf:
            a = pysam.AlignedSegment()
            a.query_name = "abc"
            a.flag = pysam.FUNMAP
            outf.write(a)

        with pysam.AlignmentFile(fn) as inf:
            read_header = inf.header

        self.compare_headers(header, read_header)
        expected_lengths = dict([(x["SN"], x["LN"]) for x in header["SQ"]])
        assert dict(zip(read_header.references, read_header.lengths)) == expected_lengths

    def test_SAM(self, tmp_path):
        self.check_read_write(tmp_path / "output.sam", "wh", self.header)

    def test_BAM(self, tmp_path):
        self.check_read_write(tmp_path / "output.bam", "wb", self.header)

    def test_CRAM(self, tmp_path):
        header = copy.copy(self.header)
        self.check_read_write(tmp_path / "output.cram", "wc", header)


class TestHeaderLargeContigs(TestHeaderWriteRead):
    """see issue 741"""

    header = {'SQ': [{'LN': 2147483647, 'SN': 'chr1'},
                     {'LN': 1584, 'SN': 'chr2'}],
              'HD': {'VN': '1.0'}}
