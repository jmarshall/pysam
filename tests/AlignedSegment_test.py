import os
import pysam
import pytest
import json
import collections
import struct
import copy
import array
from pysam import CDEL, CDIFF, CEQUAL, CINS, CMATCH, CPAD, CREF_SKIP, CSOFT_CLIP

from TestUtils import (
    dict_of_read,
    make_data_files,
    BAM_DATADIR,
    get_temp_filename,
    get_temp_context,
)


def setUpModule():
    make_data_files(BAM_DATADIR)


class ReadTest:
    def build_read(self):
        """build an example read."""

        header = pysam.AlignmentHeader.from_references(
            ["chr1", "chr2"], [10000000, 10000000]
        )

        a = pysam.AlignedSegment(header)
        a.query_name = "read_12345"
        a.query_sequence = "ATGC" * 10
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 20
        a.mapping_quality = 20
        a.cigartuples = ((0, 10), (2, 1), (0, 9), (1, 1), (0, 20))
        a.next_reference_id = 0
        a.next_reference_start = 200
        a.template_length = 167
        a.query_qualities = pysam.qualitystring_to_array("1234") * 10
        return a


class TestAlignedSegment(ReadTest):

    """tests to check if aligned read can be constructed
    and manipulated.
    """

    def check_get_aligned_pairs_combos(self, a, exp):
        def positions(exp):  return [(ppos, rpos)        for ppos, rpos, base, cigar in exp]
        def with_seq(exp):   return [(ppos, rpos, base)  for ppos, rpos, base, cigar in exp]
        def with_cigar(exp): return [(ppos, rpos, cigar) for ppos, rpos, base, cigar in exp]

        assert a.get_aligned_pairs() == positions(exp)
        assert a.get_aligned_pairs(with_seq=True) == with_seq(exp)
        assert a.get_aligned_pairs(with_cigar=True) == with_cigar(exp)
        assert a.get_aligned_pairs(with_seq=True, with_cigar=True) == exp

        exp = [(ppos, rpos, base, cigar) for ppos, rpos, base, cigar in exp if ppos is not None and rpos is not None]

        assert a.get_aligned_pairs(matches_only=True) == positions(exp)
        assert a.get_aligned_pairs(matches_only=True, with_seq=True) == with_seq(exp)
        assert a.get_aligned_pairs(matches_only=True, with_cigar=True) == with_cigar(exp)
        assert a.get_aligned_pairs(matches_only=True, with_seq=True, with_cigar=True) == exp

    def check_get_aligned_pairs_combos_without_MD(self, a, exp):
        def positions(exp): return [(ppos, rpos) for ppos, rpos, cigar in exp]

        assert a.get_aligned_pairs() == positions(exp)
        with pytest.raises(ValueError): a.get_aligned_pairs(with_seq=True)
        assert a.get_aligned_pairs(with_cigar=True) == exp
        with pytest.raises(ValueError): a.get_aligned_pairs(with_seq=True, with_cigar=True)

        exp = [(ppos, rpos, cigar) for ppos, rpos, cigar in exp if ppos is not None and rpos is not None]

        assert a.get_aligned_pairs(matches_only=True) == positions(exp)
        with pytest.raises(ValueError): a.get_aligned_pairs(matches_only=True, with_seq=True)
        assert a.get_aligned_pairs(matches_only=True, with_cigar=True) == exp
        with pytest.raises(ValueError): a.get_aligned_pairs(matches_only=True, with_seq=True, with_cigar=True)

    def testEmpty(self):

        a = pysam.AlignedSegment()
        assert a.query_name is None
        assert a.query_sequence is None
        assert pysam.qualities_to_qualitystring(a.query_qualities) is None
        assert a.flag == 0
        assert a.reference_id == -1
        assert a.mapping_quality == 0
        assert a.cigartuples is None
        assert a.get_tags() == []
        assert a.next_reference_id == -1
        assert a.next_reference_start == -1
        assert a.template_length == 0

    def testStrOfEmptyRead(self):
        a = pysam.AlignedSegment()
        s = str(a)
        assert s == "None\t0\t*\t0\t0\tNone\t*\t0\t0\tNone\tNone\t[]"

    def testSettingTagInEmptyRead(self):
        """see issue 62"""
        a = pysam.AlignedSegment()
        a.tags = (("NM", 1),)
        a.query_qualities = None
        assert a.tags == [("NM", 1),]

    def testCompare(self):
        """check comparison functions."""
        a = self.build_read()
        b = None

        assert not (a is b)
        assert not (a == b)
        assert a.compare(b) == -1

        b = self.build_read()

        assert a.compare(b) == 0
        assert b.compare(a) == 0
        assert a == b
        assert b == a
        assert not (a != b)
        assert not (b != a)

        b.tid = 1
        assert not (a == b)
        assert not (b == a)
        assert a != b
        assert b != a

    def testHashing(self):
        a = self.build_read()
        b = self.build_read()
        assert hash(a) == hash(b)
        b.tid = 1
        assert hash(a) != hash(b)

    def testUpdate(self):
        """check if updating fields affects other variable length data
        """
        a = self.build_read()
        b = self.build_read()

        # check qname
        exclude = {"query_name"}
        b.query_name = "read_123"
        assert dict_of_read(a, exclude) == dict_of_read(b, exclude)
        b.query_name = "read_12345678"
        assert dict_of_read(a, exclude) == dict_of_read(b, exclude)
        b.query_name = "read_12345"
        assert dict_of_read(a) == dict_of_read(b)

        # check cigar
        exclude = {"cigartuples"}
        b.cigartuples = ((0, 10),)
        assert dict_of_read(a, exclude) == dict_of_read(b, exclude)
        b.cigartuples = ((0, 10), (2, 1), (0, 10))
        assert dict_of_read(a, exclude) == dict_of_read(b, exclude)
        b.cigartuples = ((0, 10), (2, 1), (0, 9), (1, 1), (0, 20))
        assert dict_of_read(a) == dict_of_read(b)

        # check seq
        exclude = {"query_sequence", "query_qualities", "query_length"}
        b.query_sequence = "ATGC"
        assert dict_of_read(a, exclude) == dict_of_read(b, exclude)
        b.query_sequence = "ATGC" * 3
        assert dict_of_read(a, exclude) == dict_of_read(b, exclude)
        b.query_sequence = "ATGC" * 10
        exclude = {"query_qualities"}
        assert dict_of_read(a, exclude) == dict_of_read(b, exclude)

        # reset qual
        b = self.build_read()

        def dual(name):
            if name.endswith('is_unmapped'): return name.replace('unmapped', 'mapped')
            elif name.endswith('is_mapped'): return name.replace('mapped', 'unmapped')
            elif name.endswith('is_reverse'): return name.replace('reverse', 'forward')
            elif name.endswith('is_forward'): return name.replace('forward', 'reverse')
            else: return name

        # check flags:
        for x in (
            "is_paired",
            "is_proper_pair",
            "is_unmapped",
            "mate_is_unmapped",
            "is_reverse",
            "mate_is_reverse",
            "is_read1",
            "is_read2",
            "is_secondary",
            "is_qcfail",
            "is_duplicate",
            "is_supplementary",
        ):
            setattr(b, x, True)
            assert getattr(b, x)
            exclude = {"flag", x, dual(x)}
            assert dict_of_read(a, exclude) == dict_of_read(b, exclude)
            setattr(b, x, False)
            assert not getattr(b, x)
            assert dict_of_read(a) == dict_of_read(b)

        for x in (
            "is_mapped",
            "mate_is_mapped",
            "is_forward",
            "mate_is_forward",
        ):
            setattr(b, x, False)
            assert not getattr(b, x)
            exclude = {"flag", x, dual(x)}
            assert dict_of_read(a, exclude) == dict_of_read(b, exclude)
            setattr(b, x, True)
            assert getattr(b, x)
            assert dict_of_read(a) == dict_of_read(b)

    def testUpdate2(self):
        """issue 135: inplace update of sequence and quality score.

        This does not work as setting the sequence will erase
        the quality scores.
        """
        a = self.build_read()
        a.query_sequence = a.query_sequence[5:10]
        assert pysam.qualities_to_qualitystring(a.query_qualities) is None

        a = self.build_read()
        s = pysam.qualities_to_qualitystring(a.query_qualities)
        a.query_sequence = a.query_sequence[5:10]
        a.query_qualities = pysam.qualitystring_to_array(s[5:10])

        assert pysam.qualities_to_qualitystring(a.query_qualities) == s[5:10]

    def testClearSequence(self):
        a = pysam.AlignedSegment()
        a.query_sequence = "ATGC"
        assert a.query_sequence == "ATGC"
        a.query_sequence = None
        assert a.query_length == 0

        a.query_sequence = "ATGC"
        assert a.query_sequence == "ATGC"
        a.query_sequence = ""
        assert a.query_length == 0

        a.query_sequence = "ATGC"
        assert a.query_sequence == "ATGC"
        a.query_sequence = "*"
        assert a.query_length == 0

    def testUpdateSequenceEffects1(self):
        a = self.build_read()
        a.query_sequence = "ATGCATGC"
        a.cigarstring = "1S5M2S"
        assert a.query_alignment_sequence == "TGCAT"

        a.query_sequence = "AATTGGCC"
        assert a.query_alignment_sequence == "ATTGG"

    def testUpdateSequenceEffects2(self):
        a = self.build_read()
        a.query_sequence = "ATGCATGC"
        a.cigarstring = "1S5M2S"
        assert a.query_alignment_sequence == "TGCAT"

        a.query_sequence = "*"
        assert a.query_sequence is None
        assert a.query_alignment_sequence is None

    def testUpdateQual(self):
        """Ensure SEQ and QUAL updates leading to absent QUAL set all bytes to 0xff"""

        a = self.build_read()
        with get_temp_context("absent_qual.bam") as fname:
            with pysam.AlignmentFile(fname, "wb", header=a.header) as outf:
                a.query_sequence = "ATGC"
                outf.write(a)

                a.query_sequence = "ATGCATGCATGC"
                outf.write(a)

                a.query_sequence = "ATGCATGC"
                a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<")
                a.query_qualities = None
                outf.write(a)

            with pysam.BGZFile(fname) as f:
                # Skip BAM header
                (l_text,) = struct.unpack("<4xL", f.read(8))
                f.read(l_text)
                (n_ref,) = struct.unpack("<L", f.read(4))
                for i in range(n_ref):
                    (l_name,) = struct.unpack("<L", f.read(4))
                    f.read(l_name + 4)

                # Read each BAM record and check its qual bytes
                while True:
                    core = f.read(36)
                    if len(core) != 36: break

                    (block_size, l_read_name, n_cigar_op, l_seq) = struct.unpack("<L8xB3xH2xL12x", core)
                    data = f.read(block_size - 32)
                    qual = data[l_read_name + 4*n_cigar_op + ((l_seq+1) // 2):]

                    assert qual == b'\xff' * l_seq

    def testClearQual(self):
        a = pysam.AlignedSegment()
        a.query_sequence = "ATGC"
        a.query_qualities = pysam.qualitystring_to_array("qrst")
        a.query_qualities = None
        assert a.query_qualities is None

    def testClearQualStr(self):
        a = pysam.AlignedSegment()
        a.query_sequence = "ATGC"
        a.query_qualities_str = "qrst"
        assert a.query_qualities == pysam.qualitystring_to_array("qrst")
        assert a.query_qualities_str == "qrst"

        a.query_qualities_str = None
        assert a.query_qualities is None
        assert a.query_qualities_str is None

        a.query_qualities_str = "qrst"
        a.query_qualities_str = ""
        assert a.query_qualities is None
        assert a.query_qualities_str is None

        a.query_qualities_str = "qrst"
        a.query_qualities_str = "*"
        assert a.query_qualities is None
        assert a.query_qualities_str is None

    def testUpdateQualArrayB(self):
        a = pysam.AlignedSegment()
        a.query_sequence = "ATGC"
        a.query_qualities = array.array('B', [80, 81, 82, 83])
        assert len(a.query_qualities) == 4
        assert a.query_qualities_str == "qrst"

    def testUpdateQualArrayI(self):
        a = pysam.AlignedSegment()
        a.query_sequence = "ATGC"
        a.query_qualities = array.array('I', [80, 81, 82, 83])
        assert len(a.query_qualities) == 4
        assert a.query_qualities_str == "qrst"

    def testUpdateQualList(self):
        a = pysam.AlignedSegment()
        a.query_sequence = "ATGC"
        qual = [80, 81, 82, 83]
        a.query_qualities = qual
        qual.pop()
        assert len(a.query_qualities) == 4
        assert a.query_qualities_str == "qrst"

    def testUpdateQualString(self):
        a = pysam.AlignedSegment()
        a.query_sequence = "ATGC"
        a.query_qualities = "qrst"
        assert len(a.query_qualities) == 4
        assert a.query_qualities_str == "qrst"
        assert a.qual == "qrst"

    def testUpdateQualString2(self):
        a = pysam.AlignedSegment()
        a.query_sequence = "ATGC"
        a.query_qualities_str = "qrst"
        assert len(a.query_qualities) == 4
        assert a.query_qualities_str == "qrst"
        assert a.qual == "qrst"

    def testUpdateQualTuple(self):
        a = pysam.AlignedSegment()
        a.query_sequence = "ATGC"
        a.query_qualities = (80, 81, 82, 83)
        assert len(a.query_qualities) == 4
        assert a.query_qualities_str == "qrst"

    def testLargeRead(self):
        """build an example read."""

        a = pysam.AlignedSegment()
        a.query_name = "read_12345"
        a.query_sequence = "ATGC" * 200
        a.flag = 0
        a.reference_id = -1
        a.reference_start = 20
        a.mapping_quality = 20
        a.cigartuples = ((0, 4 * 200),)
        a.next_reference_id = 0
        a.next_reference_start = 200
        a.template_length = 167
        a.query_qualities = pysam.qualitystring_to_array("1234") * 200

        assert a

    def testUpdateTlen(self):
        """check if updating tlen works"""
        a = self.build_read()
        oldlen = a.template_length
        oldlen *= 2
        a.template_length = oldlen
        assert a.template_length == oldlen

    def testPositions(self):
        a = self.build_read()
        assert a.get_reference_positions() == \
            [
                20,
                21,
                22,
                23,
                24,
                25,
                26,
                27,
                28,
                29,
                31,
                32,
                33,
                34,
                35,
                36,
                37,
                38,
                39,
                40,
                41,
                42,
                43,
                44,
                45,
                46,
                47,
                48,
                49,
                50,
                51,
                52,
                53,
                54,
                55,
                56,
                57,
                58,
                59,
            ]

        self.check_get_aligned_pairs_combos_without_MD(
            a,
            [
                (0, 20, CMATCH),
                (1, 21, CMATCH),
                (2, 22, CMATCH),
                (3, 23, CMATCH),
                (4, 24, CMATCH),
                (5, 25, CMATCH),
                (6, 26, CMATCH),
                (7, 27, CMATCH),
                (8, 28, CMATCH),
                (9, 29, CMATCH),
                (None, 30, CDEL),
                (10, 31, CMATCH),
                (11, 32, CMATCH),
                (12, 33, CMATCH),
                (13, 34, CMATCH),
                (14, 35, CMATCH),
                (15, 36, CMATCH),
                (16, 37, CMATCH),
                (17, 38, CMATCH),
                (18, 39, CMATCH),
                (19, None, CINS),
                (20, 40, CMATCH),
                (21, 41, CMATCH),
                (22, 42, CMATCH),
                (23, 43, CMATCH),
                (24, 44, CMATCH),
                (25, 45, CMATCH),
                (26, 46, CMATCH),
                (27, 47, CMATCH),
                (28, 48, CMATCH),
                (29, 49, CMATCH),
                (30, 50, CMATCH),
                (31, 51, CMATCH),
                (32, 52, CMATCH),
                (33, 53, CMATCH),
                (34, 54, CMATCH),
                (35, 55, CMATCH),
                (36, 56, CMATCH),
                (37, 57, CMATCH),
                (38, 58, CMATCH),
                (39, 59, CMATCH),
            ],
        )

        assert a.get_reference_positions() == \
            [
                x[1]
                for x in a.get_aligned_pairs()
                if x[0] is not None and x[1] is not None
            ]
        # alen is the length of the aligned read in genome
        assert a.reference_length == a.get_aligned_pairs()[-1][0] + 1
        # aend points to one beyond last aligned base in ref
        assert a.get_reference_positions()[-1] == a.reference_end - 1

    def testFullReferencePositions(self):
        """see issue 26"""
        a = self.build_read()
        a.cigar = [(4, 30), (0, 20), (1, 3), (0, 47)]

        assert len(a.get_reference_positions(full_length=True)) == 100

    def testBlocks(self):
        a = self.build_read()
        assert a.get_blocks() == [(20, 30), (31, 40), (40, 60)]

    def test_infer_query_length(self):
        """Test infer_query_length on M|=|X|I|D|H|S cigar ops"""
        a = self.build_read()
        a.cigarstring = "40M"
        assert a.infer_query_length() == 40
        a.cigarstring = "40="
        assert a.infer_query_length() == 40
        a.cigarstring = "40X"
        assert a.infer_query_length() == 40
        a.cigarstring = "20M5I20M"
        assert a.infer_query_length() == 45
        a.cigarstring = "20M5D20M"
        assert a.infer_query_length() == 40
        a.cigarstring = "5H35M"
        assert a.infer_query_length() == 35
        a.cigarstring = "5S35M"
        assert a.infer_query_length() == 40
        a.cigarstring = "35M5H"
        assert a.infer_query_length() == 35
        a.cigarstring = "35M5S"
        assert a.infer_query_length() == 40
        a.cigarstring = None
        assert a.infer_query_length() is None

    def test_infer_read_length(self):
        """Test infer_read_length on M|=|X|I|D|H|S cigar ops"""
        a = self.build_read()
        a.cigarstring = "40M"
        assert a.infer_read_length() == 40
        a.cigarstring = "40="
        assert a.infer_read_length() == 40
        a.cigarstring = "40X"
        assert a.infer_read_length() == 40
        a.cigarstring = "20M5I20M"
        assert a.infer_read_length() == 45
        a.cigarstring = "20M5D20M"
        assert a.infer_read_length() == 40
        a.cigarstring = "5H35M"
        assert a.infer_read_length() == 40
        a.cigarstring = "5S35M"
        assert a.infer_read_length() == 40
        a.cigarstring = "35M5H"
        assert a.infer_read_length() == 40
        a.cigarstring = "35M5S"
        assert a.infer_read_length() == 40
        a.cigarstring = None
        assert a.infer_read_length() is None

    def test_get_aligned_pairs_soft_clipping(self):
        a = self.build_read()
        a.cigartuples = ((4, 2), (0, 35), (4, 3))
        self.check_get_aligned_pairs_combos_without_MD(
            a,
            [(0, None, CSOFT_CLIP), (1, None, CSOFT_CLIP)]
            + [
                (qpos, refpos, CMATCH)
                for (qpos, refpos) in zip(range(2, 2 + 35), range(20, 20 + 35))
            ]
            + [(37, None, CSOFT_CLIP), (38, None, CSOFT_CLIP), (39, None, CSOFT_CLIP)],
        )

    def test_get_aligned_pairs_hard_clipping(self):
        a = self.build_read()
        a.cigartuples = ((5, 2), (0, 35), (5, 3))
        self.check_get_aligned_pairs_combos_without_MD(
            a,
            # No seq, no seq pos
            [
                (qpos, refpos, CMATCH)
                for (qpos, refpos) in zip(range(0, 0 + 35), range(20, 20 + 35))
            ],
        )

    def test_get_aligned_pairs_skip(self):
        a = self.build_read()
        a.cigarstring = "2M100D38M"
        self.check_get_aligned_pairs_combos_without_MD(
            a,
            [(0, 20, CMATCH), (1, 21, CMATCH)]
            + [(None, refpos, CDEL) for refpos in range(22, 22 + 100)]
            + [
                (qpos, refpos, CMATCH)
                for (qpos, refpos) in zip(
                    range(2, 2 + 38), range(20 + 2 + 100, 20 + 2 + 100 + 38)
                )
            ],
        )

    def test_get_aligned_pairs_match_mismatch(self):
        a = self.build_read()
        a.cigartuples = ((7, 20), (8, 20))
        self.check_get_aligned_pairs_combos_without_MD(
            a,
            [
                (qpos, refpos, CEQUAL if qpos < 20 else CDIFF)
                for (qpos, refpos) in zip(range(0, 0 + 40), range(20, 20 + 40))
            ],
        )

    def test_get_aligned_pairs_padding(self):
        a = self.build_read()
        a.cigartuples = ((0, 1), (6, 1), (0, 1))
        # The padding operation is like an insertion into the reference.
        # See comment in test_get_aligned_pairs_padding_with_seq (below).
        self.check_get_aligned_pairs_combos_without_MD(a,
                         [(0, 20, CMATCH), (1, None, CPAD), (2, 21, CMATCH)])

    def test_get_aligned_pairs_padding_via_cigarstring(self):
        a = self.build_read()
        a.cigarstring = "1M1P1M"
        # The padding operation is like an insertion into the reference.
        # See comment in test_get_aligned_pairs_padding_with_seq (below).
        self.check_get_aligned_pairs_combos_without_MD(a,
                         [(0, 20, CMATCH), (1, None, CPAD), (2, 21, CMATCH)])

    def test_get_aligned_pairs_padding_with_seq(self):
        a = self.build_read()
        a.query_sequence = "AGT"
        a.cigarstring = "1M1P1M"
        a.set_tag("MD", "2")
        # When the reference is padded (conventionally with '*', according
        # to the SAM format specification (June 3, 2021)), as indicated by a
        # `P` CIGAR operation, this is equivalent to an insertion into the
        # reference. Thus we get the same result back as for an insertion,
        # and the reference character (if requested via with_seq=True) is
        # returned as None.
        #
        # Note that we're here assuming (as in the treatment of `N`
        # (skipped region from the reference) in build_alignment_sequence
        # in libcalignedsegment.pyx) that the MD tag will not mention the
        # region of the reference with the padding character.
        #
        # Note also that we are not dealing with a "Padded SAM" file, as
        # described in section 3.2 of the SAM format, but with the simpler
        # case (section 3.1) where a reference sequence has had '*'
        # characters inserted (by some unspecified tool) in order to make it
        # easier to specify details of an insertion using `P` in the CIGAR
        # string: "Alternatively, to describe the same alignments, we can
        # modify the reference sequence to contain pads that make room for
        # sequences inserted relative to the reference."
        self.check_get_aligned_pairs_combos(a,
            [(0, 20, "A", CMATCH), (1, None, None, CPAD), (2, 21, "T", CMATCH)])

    def test_get_aligned_pairs(self):
        a = self.build_read()
        a.query_sequence = "A" * 9
        a.cigarstring = "9M"
        a.set_tag("MD", "9")
        self.check_get_aligned_pairs_combos(
            a,
            [
                (0, 20, "A", CMATCH),
                (1, 21, "A", CMATCH),
                (2, 22, "A", CMATCH),
                (3, 23, "A", CMATCH),
                (4, 24, "A", CMATCH),
                (5, 25, "A", CMATCH),
                (6, 26, "A", CMATCH),
                (7, 27, "A", CMATCH),
                (8, 28, "A", CMATCH),
            ],
        )

        a.set_tag("MD", "4C4")
        self.check_get_aligned_pairs_combos(
            a,
            [
                (0, 20, "A", CMATCH),
                (1, 21, "A", CMATCH),
                (2, 22, "A", CMATCH),
                (3, 23, "A", CMATCH),
                (4, 24, "c", CMATCH),
                (5, 25, "A", CMATCH),
                (6, 26, "A", CMATCH),
                (7, 27, "A", CMATCH),
                (8, 28, "A", CMATCH),
            ],
        )

        a.cigarstring = "5M2D4M"
        a.set_tag("MD", "4C^TT4")
        self.check_get_aligned_pairs_combos(
            a,
            [
                (0, 20, "A", CMATCH),
                (1, 21, "A", CMATCH),
                (2, 22, "A", CMATCH),
                (3, 23, "A", CMATCH),
                (4, 24, "c", CMATCH),
                (None, 25, "T", CDEL),
                (None, 26, "T", CDEL),
                (5, 27, "A", CMATCH),
                (6, 28, "A", CMATCH),
                (7, 29, "A", CMATCH),
                (8, 30, "A", CMATCH),
            ],
        )

        a.cigarstring = "5M2D2I2M"
        a.set_tag("MD", "4C^TT2")
        self.check_get_aligned_pairs_combos(
            a,
            [
                (0, 20, "A", CMATCH),
                (1, 21, "A", CMATCH),
                (2, 22, "A", CMATCH),
                (3, 23, "A", CMATCH),
                (4, 24, "c", CMATCH),
                (None, 25, "T", CDEL),
                (None, 26, "T", CDEL),
                (5, None, None, CINS),
                (6, None, None, CINS),
                (7, 27, "A", CMATCH),
                (8, 28, "A", CMATCH),
            ],
        )

    def test_get_aligned_pairs_with_malformed_MD_tag(self):

        a = self.build_read()
        a.query_sequence = "A" * 9

        # out of range issue, see issue #560
        a.cigarstring = "64M2D85M2S"
        a.set_tag("MD", "64^TG86A0")
        with pytest.raises(AssertionError): a.get_aligned_pairs(with_seq=True)

    def test_get_aligned_pairs_skip_reference(self):
        a = self.build_read()
        a.query_sequence = "A" * 10
        a.cigarstring = "5M1N5M"
        a.set_tag("MD", "10")

        self.check_get_aligned_pairs_combos(
            a,
            [
                (0, 20, "A", CMATCH),
                (1, 21, "A", CMATCH),
                (2, 22, "A", CMATCH),
                (3, 23, "A", CMATCH),
                (4, 24, "A", CMATCH),
                (None, 25, None, CREF_SKIP),
                (5, 26, "A", CMATCH),
                (6, 27, "A", CMATCH),
                (7, 28, "A", CMATCH),
                (8, 29, "A", CMATCH),
                (9, 30, "A", CMATCH),
            ]
        )

    def test_equivalence_matches_only_and_with_seq(self):
        a = self.build_read()
        a.query_sequence = "ACGT" * 2
        a.cigarstring = "4M1D4M"
        a.set_tag("MD", "4^x4")
        self.check_get_aligned_pairs_combos(
            a,
            list(zip(range(0, 4), range(20, 24), "ACGT", [CMATCH] * 4))
            + [(None, 24, "x", CDEL)]
            + list(zip(range(4, 8), range(25, 29), "ACGT", [CMATCH] * 4)),
        )

        a = self.build_read()
        a.query_sequence = "ACGT" * 2
        a.cigarstring = "4M1N4M"
        a.set_tag("MD", "8")
        self.check_get_aligned_pairs_combos(
            a,
            list(zip(range(0, 4), range(20, 24), "ACGT", [CMATCH] * 4))
            + [(None, 24, None, 3)]
            + list(zip(range(4, 8), range(25, 29), "ACGT", [CMATCH] * 4)),
        )

    def test_get_aligned_pairs_lowercase_md(self):
        a = self.build_read()
        a.query_sequence = "A" * 10
        a.cigarstring = "10M"
        a.set_tag("MD", "5g4")
        assert a.get_aligned_pairs(with_seq=True) == \
            [
                (0, 20, "A"),
                (1, 21, "A"),
                (2, 22, "A"),
                (3, 23, "A"),
                (4, 24, "A"),
                (5, 25, "g"),
                (6, 26, "A"),
                (7, 27, "A"),
                (8, 28, "A"),
                (9, 29, "A"),
            ]

    def test_get_aligned_pairs_uppercase_md(self):
        a = self.build_read()
        a.query_sequence = "A" * 10
        a.cigarstring = "10M"
        a.set_tag("MD", "5G4")
        assert a.get_aligned_pairs(with_seq=True) == \
            [
                (0, 20, "A"),
                (1, 21, "A"),
                (2, 22, "A"),
                (3, 23, "A"),
                (4, 24, "A"),
                (5, 25, "g"),
                (6, 26, "A"),
                (7, 27, "A"),
                (8, 28, "A"),
                (9, 29, "A"),
            ]

    def test_get_aligned_pairs_1character_md(self):
        a = self.build_read()
        a.query_sequence = "A" * 7
        a.cigarstring = "7M"
        a.set_tag("MD", "7", value_type="A")
        assert a.get_aligned_pairs(with_seq=True) == \
            [
                (0, 20, "A"),
                (1, 21, "A"),
                (2, 22, "A"),
                (3, 23, "A"),
                (4, 24, "A"),
                (5, 25, "A"),
                (6, 26, "A"),
            ]

    def test_get_aligned_pairs_bad_type_md(self):
        a = self.build_read()
        a.query_sequence = "A" * 7
        a.cigarstring = "7M"
        a.set_tag("MD", 7)
        with pytest.raises(TypeError):
            a.get_aligned_pairs(with_seq=True)

    def testNoSequence(self):
        """issue 176: retrieving length without query sequence
        with soft-clipping.
        """
        a = self.build_read()
        a.query_sequence = None
        a.cigarstring = "20M"
        assert a.query_alignment_length == 20
        a.cigarstring = "20M1S"
        assert a.query_alignment_length == 20
        a.cigarstring = "20M1H"
        assert a.query_alignment_length == 20
        a.cigarstring = "1S20M"
        assert a.query_alignment_length == 20
        a.cigarstring = "1H20M"
        assert a.query_alignment_length == 20
        a.cigarstring = "1S20M1S"
        assert a.query_alignment_length == 20
        a.cigarstring = "1H20M1H"
        assert a.query_alignment_length == 20

    def test_query_length_is_limited(self):
        a = self.build_read()
        a.query_name = "A" * 1
        a.query_name = "A" * 254
        with pytest.raises(ValueError): a.query_name = "A" * 255

    def test_header_accessible(self):
        a = self.build_read()
        assert isinstance(a.header, pysam.AlignmentHeader)

    def test_bin_values_for_unmapped_reads_ignore_length(self):
        a = self.build_read()
        # use a long read
        a.cigarstring = "2000000M"
        assert a.bin == 9
        # changing unmapped flag changes bin because length is 0
        a.is_unmapped = True
        assert a.is_unmapped
        assert not a.is_mapped
        assert a.bin == 4681

        # unmapped read without chromosomal location
        a.reference_start = -1
        assert a.reference_start == -1
        assert a.bin == 4680

    def test_bin_values_for_mapped_reads_are_updated(self):
        a = self.build_read()
        a.pos = 20000
        assert not a.is_unmapped
        assert a.is_mapped
        assert a.bin == 4682

        # updating length updates bin
        a.cigarstring = "2000000M"
        assert a.bin == 9

        # updating length updates bin
        a.cigarstring = "20M"
        assert a.bin == 4682

        # updating length updates bin
        a.reference_start = 2000000
        assert a.bin == 4803


class TestTidMapping(ReadTest):
    def test_reference_name_can_be_set_to_none(self):
        a = self.build_read()
        a.reference_name = None
        assert a.reference_name is None
        assert a.reference_id == -1

    def test_reference_name_can_be_set_to_asterisk(self):
        a = self.build_read()
        a.reference_name = "*"
        assert a.reference_name is None
        assert a.reference_id == -1

    def test_reference_name_can_be_set_to_chromosome(self):
        a = self.build_read()
        a.reference_name = "chr1"
        assert a.reference_name == "chr1"
        assert a.reference_id == 0

    def test_reference_name_can_not_be_set_to_unknown_chromosome(self):
        a = self.build_read()
        with pytest.raises(ValueError): a.reference_name = "chrX"

    def test_tid_can_be_set_to_missing(self):
        a = self.build_read()
        a.reference_id = -1
        assert a.reference_id == -1
        assert a.reference_name is None

    def test_tid_can_be_set_to_missing_without_header(self):
        a = pysam.AlignedSegment()
        a.reference_id = -1
        assert a.reference_id == -1
        assert a.reference_name is None

    def test_tid_can_be_set_without_header(self):
        a = pysam.AlignedSegment()
        a.reference_id = 1
        with pytest.raises(ValueError): a.reference_name

    def test_tid_can_be_set_to_chromosome(self):
        a = self.build_read()
        a.reference_id = 0
        assert a.reference_id == 0
        assert a.reference_name == "chr1"

    def test_tid_can_not_be_set_to_unknown_chromosome(self):
        a = self.build_read()
        with pytest.raises(ValueError): a.reference_id = 2

    def test_unmapped_tid_is_asterisk_in_output(self):
        a = self.build_read()
        a.reference_id = -1
        assert a.to_string().split("\t")[2] == "*"


class TestNextTidMapping(ReadTest):
    def test_next_reference_name_can_be_set_to_none(self):
        a = self.build_read()
        a.next_reference_name = None
        assert a.next_reference_name is None
        assert a.next_reference_id == -1

    def test_next_reference_name_can_be_set_to_asterisk(self):
        a = self.build_read()
        a.next_reference_name = "*"
        assert a.next_reference_name is None
        assert a.next_reference_id == -1

    def test_next_reference_name_can_be_set_to_chromosome(self):
        a = self.build_read()
        a.next_reference_name = "chr1"
        assert a.next_reference_name == "chr1"
        assert a.next_reference_id == 0

    def test_next_reference_name_can_not_be_set_to_unknown_chromosome(self):
        a = self.build_read()
        with pytest.raises(ValueError): a.next_reference_name = "chrX"

    def test_next_tid_can_be_set_to_missing(self):
        a = self.build_read()
        a.next_reference_id = -1
        assert a.next_reference_id == -1
        assert a.next_reference_name is None

    def test_next_tid_can_be_set_to_equal(self):
        a = self.build_read()
        a.reference_name = "chr1"
        a.next_reference_name = "="
        assert a.next_reference_id == a.reference_id
        assert a.next_reference_name == a.reference_name
        assert a.to_string().split("\t")[6] == "="

    def test_next_tid_can_be_set_to_missing_without_header(self):
        a = pysam.AlignedSegment()
        a.next_reference_id = -1
        assert a.next_reference_id == -1
        assert a.next_reference_name is None

    def test_next_tid_can_be_set_without_header(self):
        a = pysam.AlignedSegment()
        a.next_reference_id = 1
        with pytest.raises(ValueError):  a.next_reference_name

    def test_next_tid_can_be_set_to_chromosome(self):
        a = self.build_read()
        a.next_reference_id = 0
        assert a.next_reference_id == 0
        assert a.next_reference_name == "chr1"

    def test_next_tid_can_not_be_set_to_unknown_chromosome(self):
        a = self.build_read()
        with pytest.raises(ValueError): a.next_reference_id = 2

    def test_next_unmapped_tid_is_asterisk_in_output(self):
        a = self.build_read()
        a.next_reference_id = -1
        assert a.to_string().split("\t")[6] == "*"


class TestCigar(ReadTest):
    def testCigarString(self):
        r = self.build_read()
        assert r.cigarstring == "10M1D9M1I20M"
        r.cigarstring = "20M10D20M"
        assert r.cigartuples == [(0, 20), (2, 10), (0, 20)]
        # unsetting cigar string
        r.cigarstring = None
        assert r.cigarstring is None

        r.cigarstring = "40M"
        assert r.cigartuples == [(0, 40)]
        r.cigarstring = ""
        assert r.cigarstring is None

        r.cigarstring = "40M"
        assert r.cigartuples == [(0, 40)]
        r.cigarstring = "*"
        assert r.cigarstring is None

    def testCigar(self):
        r = self.build_read()
        assert r.cigartuples == [(0, 10), (2, 1), (0, 9), (1, 1), (0, 20)]
        # unsetting cigar string
        r.cigartuples = None
        assert r.cigartuples is None


class TestCigarStats(ReadTest):
    def testStats(self):
        a = self.build_read()

        a.cigarstring = None
        assert [list(x) for x in a.get_cigar_stats()] == \
            [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

        a.cigarstring = "10M"
        assert [list(x) for x in a.get_cigar_stats()] == \
            [[10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

        a.cigarstring = "10M2I2M"
        assert [list(x) for x in a.get_cigar_stats()] == \
            [[12, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0], [2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

        for i, x in enumerate("MIDNSHP=X"):
            a.cigarstring = "2{}".format(x)
            expected = [[0] * 11, [0] * 11]
            expected[0][i] = 2
            expected[1][i] = 1
            assert [list(x) for x in a.get_cigar_stats()] == expected

        for i in range(1, 100):
            cigarstring = "".join("10{}".format(x)
                                  for x in iter("MIDNSHP=X")) * i
            a.cigarstring = cigarstring
            assert a.cigarstring == cigarstring
            expected = [[i * 10 for j in range(len("MIDNSHP=X"))] + [0, 0],
                        [i for j in range(len("MIDNSHP=X"))] + [0, 0]]
            obtained = [list(x) for x in a.get_cigar_stats()]
            assert obtained == expected

        a.cigarstring = "10M"
        a.set_tag("NM", 5)
        assert [list(x) for x in a.get_cigar_stats()] == \
            [[10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

        a.cigarstring = None
        assert [list(x) for x in a.get_cigar_stats()] == \
            [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]


class TestAlignedPairs:
    filename = os.path.join(BAM_DATADIR, "example_aligned_pairs.bam")

    def testReferenceBases(self):
        """reference bases should always be the same nucleotide
        """
        reference_bases = collections.defaultdict(list)
        with pysam.AlignmentFile(self.filename) as inf:
            for c in inf.pileup():
                for r in c.pileups:
                    for read, ref, base in r.alignment.get_aligned_pairs(with_seq=True):
                        if ref is None:
                            continue
                        reference_bases[ref].append(base.upper())

        for x, y in reference_bases.items():
            assert len(set(y)) == 1


class TestBaseModifications:
    def testChebi(self):
        """reference bases should always be the same nucleotide
        """
        filename = os.path.join(BAM_DATADIR, "MM-chebi.bam")
        expect = {
            ("C", 0, "m"): [(6, 102), (17, 128), (20, 153), (31, 179), (34, 204)],
            ("N", 0, "n"): [(15, 212)],
            ("C", 0, 76792): [(19, 161), (34, 33)],
        }

        with pysam.AlignmentFile(filename, check_sq=False) as inf:
            r = next(iter(inf))
            assert r.modified_bases == expect

    def testDouble(self):
        """reference bases should always be the same nucleotide
        """
        filename = os.path.join(BAM_DATADIR, "MM-double.bam")
        expect = {
            ("G", 1, "m"): [(1, 115), (12, 141), (13, 166), (22, 192)],
            ("C", 0, "m"): [(7, 128), (30, 153), (31, 179)],
            ("G", 0, "o"): [(13, 102)],
        }

        with pysam.AlignmentFile(filename, check_sq=False) as inf:
            r = next(iter(inf))
            assert r.modified_bases == expect

    def testExplicit(self):
        """reference bases should always be the same nucleotide
        """
        filename = os.path.join(BAM_DATADIR, "MM-explicit.bam")
        expected_output = [
            {("C", 0, "m"): [(9, 200), (10, 50), (14, 160)], ("C", 0, "h"): [(9, 10), (10, 170), (14, 20)]},
            {("C", 0, "m"): [(9, 200), (10, 50), (13, 10), (14, 160), (16, 10)],
             ("C", 0, "h"): [(9, 10), (10, 170), (13, 5), (14, 20), (16, 5)]},
            {("C", 0, "m"): [(9, 200), (14, 160)], ("C", 0, "h"): [(9, 10), (10, 170), (13, 5), (14, 20), (16, 5)]},
        ]

        with pysam.AlignmentFile(filename, check_sq=False) as inf:
            for r, expected in zip(inf, expected_output):
                assert r.modified_bases == expected

    def testMulti(self):
        """reference bases should always be the same nucleotide
        """
        filename = os.path.join(BAM_DATADIR, "MM-multi.bam")
        expect = {
            "r1": {
                ("C", 0, "m"): [(6, 128), (17, 153), (20, 179), (31, 204), (34, 230)],
                ("N", 0, "n"): [(15, 215), (18, 240)],
                ("C", 0, "h"): [(19, 159), (34, 6)],
            },
            "r2": {
                ("C", 0, "m"): [
                    (6, 77),
                    (17, 103),
                    (19, 128),
                    (20, 154),
                    (31, 179),
                    (34, 204),
                ],
                ("C", 0, "h"): [
                    (6, 159),
                    (17, 133),
                    (19, 108),
                    (20, 82),
                    (31, 57),
                    (34, 31),
                ],
                ("N", 0, "n"): [(15, 240)],
            },
        }

        with pysam.AlignmentFile(filename, check_sq=False) as inf:
            for r in inf:
                assert r.modified_bases == expect[r.query_name]

    def testOrient(self):
        """reference bases should always be the same nucleotide
        """
        filename = os.path.join(BAM_DATADIR, "MM-orient.bam")
        expect = {
            "top-fwd": [
                {("C", 0, "m"): [(7, 128), (30, 153), (31, 179)]},
                {("C", 0, "m"): [(7, 128), (30, 153), (31, 179)]},
            ],
            "top-rev": [
                {("C", 1, "m"): [(4, 179), (5, 153), (28, 128)]},
                {("C", 0, "m"): [(31, 179), (30, 153), (7, 128)]},
            ],
            "bot-fwd": [
                {("G", 1, "m"): [(1, 115), (2, 141), (18, 166), (23, 192)]},
                {("G", 1, "m"): [(1, 115), (2, 141), (18, 166), (23, 192)]},
            ],
            "bot-rev": [
                {("G", 0, "m"): [(12, 192), (17, 166), (33, 141), (34, 115)]},
                {("G", 1, "m"): [(23, 192), (18, 166), (2, 141), (1, 115)]},
            ],
        }

        with pysam.AlignmentFile(filename, check_sq=False) as inf:
            for r in inf:
                assert r.modified_bases == expect[r.query_name][0]
                assert r.modified_bases_forward == expect[r.query_name][1]
                for (B, s, _), mods in r.modified_bases.items():
                    C = B.translate(str.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx"))
                    for pos, _ in mods:
                        if r.is_reverse:
                            if s == 1:
                                assert C == r.query_sequence[pos], r.to_string()
                            else:
                                assert C == r.query_sequence[pos], r.to_string()
                        else:
                            if s == 0:
                                assert B == r.query_sequence[pos], r.to_string()
                            else:
                                assert B == r.query_sequence[pos], r.to_string()


class TestTags(ReadTest):
    def testMissingTag(self):
        a = self.build_read()
        with pytest.raises(KeyError): a.get_tag("XP")

    def testEmptyTag(self):
        a = self.build_read()
        with pytest.raises(KeyError): a.get_tag("XT")

    def testSetTag(self):
        a = self.build_read()
        assert not a.has_tag("NM")
        a.set_tag("NM", 2)
        assert a.has_tag("NM")
        assert a.get_tag("NM") == 2
        a.set_tag("NM", 3)
        assert a.get_tag("NM") == 3
        a.set_tag("NM", None)
        assert not a.has_tag("NM")
        # check if deleting a non-existing tag is fine
        a.set_tag("NM", None)
        a.set_tag("NM", None)

    def testArrayTags(self):
        read = self.build_read()
        supported_dtypes = "bhBHf"
        unsupported_dtypes = "lLd"

        for dtype in supported_dtypes:
            key = "F" + dtype
            read.set_tag(key, array.array(dtype, range(10)))
            ary = read.get_tag(key)

        for dtype in unsupported_dtypes:
            key = "F" + dtype
            with pytest.raises(ValueError): read.set_tag(key, array.array(dtype, range(10)))

    def testAddTagsType(self):
        a = self.build_read()
        a.tags = None
        assert a.tags == []

        a.setTag("X1", 5.0)
        a.setTag("X2", "5.0")
        a.setTag("X3", 5)

        assert sorted(a.tags) == sorted([("X1", 5.0), ("X2", "5.0"), ("X3", 5)])

        # test setting float for int value
        a.setTag("X4", 5, value_type="d")
        assert sorted(a.tags) == sorted([("X1", 5.0), ("X2", "5.0"), ("X3", 5), ("X4", 5.0)])

        # test setting int for float value - the
        # value will be rounded.
        a.setTag("X5", 5.2, value_type="i")
        assert sorted(a.tags) == sorted([("X1", 5.0), ("X2", "5.0"), ("X3", 5), ("X4", 5.0), ("X5", 5)])

        # test setting invalid type code
        with pytest.raises(ValueError): a.set_tag("X6", 5.2, "g")

    def testTagsUpdatingFloat(self):
        a = self.build_read()
        a.tags = [("NM", 1), ("RG", "L1"), ("PG", "P1"), ("XT", "U")]

        assert a.tags == [("NM", 1), ("RG", "L1"), ("PG", "P1"), ("XT", "U")]
        a.tags += [("XC", 5.0)]
        assert a.tags == [("NM", 1), ("RG", "L1"), ("PG", "P1"), ("XT", "U"), ("XC", 5.0)]

    def testAddTags(self):
        a = self.build_read()
        a.tags = [("NM", 1), ("RG", "L1"), ("PG", "P1"), ("XT", "U")]

        assert sorted(a.tags) == sorted([("NM", 1), ("RG", "L1"), ("PG", "P1"), ("XT", "U")])

        a.setTag("X1", "C")
        assert sorted(a.tags) == sorted([("X1", "C"), ("NM", 1), ("RG", "L1"), ("PG", "P1"), ("XT", "U"),])
        a.setTag("X2", 5)
        assert sorted(a.tags) == sorted(
                [
                    ("X2", 5),
                    ("X1", "C"),
                    ("NM", 1),
                    ("RG", "L1"),
                    ("PG", "P1"),
                    ("XT", "U"),
                ]
            )
        # add with replacement
        a.setTag("X2", 10)
        assert sorted(a.tags) == sorted(
                [
                    ("X2", 10),
                    ("X1", "C"),
                    ("NM", 1),
                    ("RG", "L1"),
                    ("PG", "P1"),
                    ("XT", "U"),
                ]
            )

        # add without replacement
        a.setTag("X2", 5, replace=False)
        assert sorted(a.tags) == sorted(
                [
                    ("X2", 10),
                    ("X1", "C"),
                    ("X2", 5),
                    ("NM", 1),
                    ("RG", "L1"),
                    ("PG", "P1"),
                    ("XT", "U"),
                ]
            )

    def testTagParsing(self):
        """test for tag parsing

        see http://groups.google.com/group/pysam-user-group/browse_thread/thread/67ca204059ea465a
        """
        samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex8.bam"), "rb")

        for entry in samfile:
            before = entry.get_tags()
            entry.set_tags(before)
            after = entry.get_tags()
            assert after == before

    def testMDTagMissing(self):
        a = self.build_read()
        with pytest.raises(ValueError): a.get_reference_sequence()

    def testMDTagMissingCigar(self):
        a = self.build_read()
        a.set_tag("MD", "5")
        a.cigartuples = None
        with pytest.raises(ValueError): a.get_reference_sequence()

    def testMDTagMissingSeq(self):
        a = self.build_read()
        a.set_tag("MD", "5")
        a.query_sequence = None
        with pytest.raises(ValueError): a.get_reference_sequence()

    def testMDTagMatchOnly(self):
        a = self.build_read()

        # Substitutions only
        a.cigarstring = "21M"
        a.query_sequence = "A" * 21
        a.set_tag("MD", "5C0T0G05C0G0T5")
        assert a.get_reference_sequence() == "AAAAActgAAAAAcgtAAAAA"

        a.cigarstring = "21M"
        a.query_sequence = "A" * 21
        a.set_tag("MD", "5CTG5CGT5")
        assert a.get_reference_sequence() == "AAAAActgAAAAAcgtAAAAA"

        a.cigarstring = "11M"
        a.query_sequence = "A" * 11
        a.set_tag("MD", "CTG5CGT")
        assert a.get_reference_sequence() == "ctgAAAAAcgt"

    def testMDTagInsertions(self):
        a = self.build_read()

        # insertions are silent in the reference sequence
        a.cigarstring = "5M1I5M"
        a.query_sequence = "A" * 5 + "C" + "A" * 5
        a.set_tag("MD", "10")
        assert a.get_reference_sequence() == "A" * 10

        a.cigarstring = "1I10M"
        a.query_sequence = "C" * 1 + "A" * 10
        assert a.get_reference_sequence() == "A" * 10

        a.cigarstring = "10M1I"
        a.query_sequence = "A" * 10 + "C" * 1
        assert a.get_reference_sequence() == "A" * 10

    def testMDTagDeletions(self):
        a = self.build_read()

        a.cigarstring = "5M1D5M"
        a.query_sequence = "A" * 10
        a.set_tag("MD", "5^C5")
        assert a.get_reference_sequence() == "A" * 5 + "C" + "A" * 5

        a.cigarstring = "5M3D5M"
        a.query_sequence = "A" * 10
        a.set_tag("MD", "5^CCC5")
        assert a.get_reference_sequence() == "A" * 5 + "C" * 3 + "A" * 5

    def testMDTagRefSkipping(self):
        a = self.build_read()

        a.cigarstring = "5M1N5M"
        a.query_sequence = "A" * 10
        a.set_tag("MD", "10")
        assert a.get_reference_sequence() == "A" * 10

        a.cigarstring = "5M3N5M"
        a.query_sequence = "A" * 10
        a.set_tag("MD", "10")
        assert a.get_reference_sequence() == "A" * 10

    def testMDTagSoftClipping(self):
        a = self.build_read()

        # softclipping
        a.cigarstring = "5S5M1D5M5S"
        a.query_sequence = "G" * 5 + "A" * 10 + "G" * 5
        a.set_tag("MD", "5^C5")
        assert a.get_reference_sequence() == "A" * 5 + "C" + "A" * 5

        # all together
        a.cigarstring = "5S5M1D5M1I5M5S"
        a.query_sequence = "G" * 5 + "A" * 16 + "G" * 5
        a.set_tag("MD", "2C2^T10")
        assert a.get_reference_sequence() == "AAcAATAAAAAAAAAA"

    def testMDTagComplex(self):
        a = self.build_read()

        a.cigarstring = "5S5M1I2D5M5S"
        a.query_sequence = "G" * 5 + "A" * 11 + "G" * 5
        a.set_tag("MD", "2C2^TC5")
        assert a.get_reference_sequence() == "AAcAATCAAAAA"

        a.cigarstring = "5S5M2D1I5M5S"
        a.query_sequence = "G" * 5 + "A" * 11 + "G" * 5
        a.set_tag("MD", "2C2^TC5")
        assert a.get_reference_sequence() == "AAcAATCAAAAA"

        # insertion in reference overlapping deletion in reference
        # read: AACCCCA---AAA
        # ref:  AA----AGGGAAA
        a.cigarstring = "2M4I1M3D3M"
        a.set_tag("MD", "3^GGG3")
        a.query_sequence = "AACCCCAAAA"
        assert a.get_reference_sequence() == "AAAGGGAAA"

        a.cigarstring = "5M2D2I2M"
        a.set_tag("MD", "4C^TT2")
        a.query_sequence = "A" * 9
        assert a.get_reference_sequence() == "AAAAcTTAA"

    def testArrayTagValues(self):

        r = self.build_read()

        def c(r, l):
            r.tags = [("ZM", l)]
            assert list(r.opt("ZM")) == list(l)

        # signed integers
        c(r, (-1, 1))
        c(r, (-1, 100))
        c(r, (-1, 200))
        c(r, (-1, 1000))
        c(r, (-1, 30000))
        c(r, (-1, 50000))
        c(r, (1, -1))
        c(r, (1, -100))
        c(r, (1, -200))
        c(r, (1, -1000))
        c(r, (1, -30000))
        c(r, (1, -50000))

        # unsigned integers
        c(r, (1, 100))
        c(r, (1, 1000))
        c(r, (1, 10000))
        c(r, (1, 100000))

        # floats
        c(r, (1.0, 100.0))

    def testLongTags(self):
        """see issue 115"""

        r = self.build_read()
        rg = "HS2000-899_199.L3"
        tags = [
            ("XC", 85),
            ("XT", "M"),
            ("NM", 5),
            ("SM", 29),
            ("AM", 29),
            ("XM", 1),
            ("XO", 1),
            ("XG", 4),
            ("MD", "37^ACCC29T18"),
            (
                "XA",
                "5,+11707,36M1I48M,2;21,-48119779,46M1I38M,2;hs37d5,-10060835,40M1D45M,3;5,+11508,36M1I48M,3;hs37d5,+6743812,36M1I48M,3;19,-59118894,46M1I38M,3;4,-191044002,6M1I78M,3;",
            ),
        ]  # noqa

        r.tags = tags
        r.tags += [("RG", rg)] * 100
        tags += [("RG", rg)] * 100

        assert r.tags == tags

    def testNegativeIntegers(self):
        x = -2
        aligned_read = self.build_read()
        aligned_read.tags = [("XD", int(x))]
        assert aligned_read.opt("XD") == x
        # print (aligned_read.tags)

    def testNegativeIntegersWrittenToFile(self):
        r = self.build_read()
        x = -2
        r.tags = [("XD", x)]
        with get_temp_context("negative_integers.bam") as fn:
            with pysam.AlignmentFile(
                fn, "wb", referencenames=("chr1",), referencelengths=(1000,)
            ) as outf:
                outf.write(r)
            with pysam.AlignmentFile(fn) as inf:
                r = next(inf)
            assert r.tags == [("XD", x)]


class TestCopy(ReadTest):
    def testCopy(self):
        a = self.build_read()
        b = copy.copy(a)
        # check if a and be are the same
        assert a == b

        # check if they map to different objects
        a.query_name = "ReadA"
        b.query_name = "ReadB"
        assert a.query_name == "ReadA"
        assert b.query_name == "ReadB"

    def testDeepCopy(self):
        a = self.build_read()
        b = copy.deepcopy(a)
        # check if a and be are the same
        assert a == b

        # check if they map to different objects
        a.query_name = "ReadA"
        b.query_name = "ReadB"
        assert a.query_name == "ReadA"
        assert b.query_name == "ReadB"


class TestSetTagGetTag(ReadTest):
    def check_tag(self, tag, value, value_type, alt_value_type=None):
        a = self.build_read()
        a.set_tag(tag, value, value_type=value_type)
        v, t = a.get_tag(tag, with_value_type=True)
        assert v == value

        if alt_value_type:
            assert t == alt_value_type
        else:
            assert t == value_type

    def test_set_tag_with_A(self):
        self.check_tag("TT", "x", value_type="A")

    def test_set_tag_with_a(self):
        self.check_tag("TT", "x", value_type="a", alt_value_type="A")

    def test_set_tag_with_C(self):
        self.check_tag("TT", 12, value_type="C")

    def test_set_tag_with_c(self):
        self.check_tag("TT", 12, value_type="c")

    def test_set_tag_with_S(self):
        self.check_tag("TT", 12, value_type="S")

    def test_set_tag_with_s(self):
        self.check_tag("TT", 12, value_type="s")

    def test_set_tag_with_I(self):
        self.check_tag("TT", 12, value_type="I")

    def test_set_tag_with_i(self):
        self.check_tag("TT", 12, value_type="i")

    def test_set_tag_with_f(self):
        self.check_tag("TT", 2.5, value_type="f")

    def test_set_tag_with_d(self):
        self.check_tag("TT", 2.5, value_type="d")

    def test_set_tag_with_H(self):
        self.check_tag("TT", "AE12", value_type="H")

    def test_set_tag_with_automated_type_detection(self):
        self.check_tag("TT", -(1 << 7), value_type=None, alt_value_type="c")
        self.check_tag("TT", -(1 << 7) - 1, value_type=None, alt_value_type="s")
        self.check_tag("TT", -(1 << 15), value_type=None, alt_value_type="s")
        self.check_tag("TT", -(1 << 15) - 1, value_type=None, alt_value_type="i")
        self.check_tag("TT", -(1 << 31), value_type=None, alt_value_type="i")
        with pytest.raises(ValueError):
            self.check_tag("TT", -(1 << 31) - 1, value_type=None, alt_value_type="i")

        self.check_tag("TT", (1 << 8) - 1, value_type=None, alt_value_type="C")
        self.check_tag("TT", (1 << 8), value_type=None, alt_value_type="S")
        self.check_tag("TT", (1 << 16) - 1, value_type=None, alt_value_type="S")
        self.check_tag("TT", (1 << 16), value_type=None, alt_value_type="I")
        self.check_tag("TT", (1 << 32) - 1, value_type=None, alt_value_type="I")
        with pytest.raises(ValueError):
            self.check_tag("TT", (1 << 32), value_type=None, alt_value_type="I")

    def test_set_tag_invalid_value_type(self):
        with pytest.raises(ValueError):
            self.check_tag("TT", "abc", value_type="#")

    def test_set_array_tag_invalid_value_type(self):
        with pytest.raises(ValueError):
            self.check_tag("TT", array.array('I', range(4)), value_type='#')

    def test_set_array_tag_invalid_typecode(self):
        with pytest.raises(ValueError):
            self.check_tag("TT", array.array('L', range(4)), value_type=None)


class TestSetTagsGetTag(TestSetTagGetTag):
    def check_tag(self, tag, value, value_type, alt_value_type=None):
        a = self.build_read()
        a.set_tags([(tag, value, value_type)])
        v, t = a.get_tag(tag, with_value_type=True)
        if alt_value_type:
            assert t == alt_value_type
        else:
            assert t == value_type
        assert v == value


def test_cigar_enums_are_defined():
    assert pysam.CIGAR_OPS.CMATCH == 0
    assert pysam.CIGAR_OPS.CINS == 1
    assert pysam.CIGAR_OPS.CDEL == 2
    assert pysam.CIGAR_OPS.CREF_SKIP == 3
    assert pysam.CIGAR_OPS.CSOFT_CLIP == 4
    assert pysam.CIGAR_OPS.CHARD_CLIP == 5
    assert pysam.CIGAR_OPS.CPAD == 6
    assert pysam.CIGAR_OPS.CEQUAL == 7
    assert pysam.CIGAR_OPS.CDIFF == 8
    assert pysam.CIGAR_OPS.CBACK == 9


def test_sam_flags_are_defined():
    assert pysam.SAM_FLAGS.FPAIRED == 1
    assert pysam.SAM_FLAGS.FPROPER_PAIR == 2
    assert pysam.SAM_FLAGS.FUNMAP == 4
    assert pysam.SAM_FLAGS.FMUNMAP == 8
    assert pysam.SAM_FLAGS.FREVERSE == 16
    assert pysam.SAM_FLAGS.FMREVERSE == 32
    assert pysam.SAM_FLAGS.FREAD1 == 64
    assert pysam.SAM_FLAGS.FREAD2 == 128
    assert pysam.SAM_FLAGS.FSECONDARY == 256
    assert pysam.SAM_FLAGS.FQCFAIL == 512
    assert pysam.SAM_FLAGS.FDUP == 1024
    assert pysam.SAM_FLAGS.FSUPPLEMENTARY == 2048


class TestBuildingReadsWithoutHeader:
    def build_read(self):
        """build an example read, but without header information."""

        a = pysam.AlignedSegment()
        a.query_name = "read_12345"
        a.query_sequence = "ATGC" * 10
        a.flag = 0
        a.reference_id = -1
        a.reference_start = 20
        a.mapping_quality = 20
        a.cigartuples = ((0, 10), (2, 1), (0, 9), (1, 1), (0, 20))
        a.next_reference_id = 0
        a.next_reference_start = 200
        a.template_length = 167
        a.query_qualities = pysam.qualitystring_to_array("1234") * 10
        # todo: create tags
        return a

    def test_read_can_be_constructed_without_header(self):
        read = self.build_read()
        assert read.query_name == "read_12345"

    def test_reference_id_can_be_set(self):
        read = self.build_read()
        read.reference_id = 2
        assert read.reference_id == 2

    def test_reference_name_is_not_available(self):
        read = self.build_read()
        with pytest.raises(ValueError): read.reference_name = "chr2"

    def test_read_can_be_written_to_file(self):
        tmpfilename = get_temp_filename(".bam")
        with pysam.AlignmentFile(
            tmpfilename,
            "wb",
            reference_names=["chr1", "chr2", "chr3"],
            reference_lengths=[1000, 2000, 3000],
        ) as outf:
            read = self.build_read()
            read.reference_id = 2
            outf.write(read)

        stdout = pysam.samtools.view(tmpfilename)
        chromosome = stdout.split("\t")[2]
        assert chromosome == "chr3"
        os.unlink(tmpfilename)


class TestForwardStrandValues(ReadTest):
    def test_sequence_is_complemented(self):
        a = self.build_read()
        a.is_reverse = False
        fwd_seq = a.query_sequence

        rev_seq = fwd_seq.translate(str.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx"))[::-1]
        assert a.get_forward_sequence() == fwd_seq
        a.is_reverse = True
        assert a.query_sequence == fwd_seq
        assert a.get_forward_sequence() == rev_seq

    def test_ambiguous_bases_are_complemented(self):
        a = self.build_read()
        ambiguity = "ABCDGHKMNRSTVWY"
        reviguity = "TVGHCDMKNYSABWR"[::-1]

        a.query_sequence = ambiguity
        a.is_reverse = False
        assert a.get_forward_sequence() == ambiguity
        a.is_reverse = True
        assert a.get_forward_sequence() == reviguity

    def test_qualities_are_complemented(self):
        a = self.build_read()
        a.is_reverse = False
        fwd_qual = a.query_qualities
        rev_qual = fwd_qual[::-1]
        assert a.get_forward_qualities() == fwd_qual
        a.is_reverse = True
        assert a.query_qualities == fwd_qual
        assert a.get_forward_qualities() == rev_qual


class TestExportImport(ReadTest):
    def test_string_export(self):
        a = self.build_read()
        assert a.to_string() == (
            "read_12345\t0\tchr1\t21\t20\t10M1D9M1I20M\t=\t201\t167\t"
            "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\t1234123412341234123412341234123412341234"
        )

    def test_string_export_import_str_without_tags(self):
        a = self.build_read()
        a.tags = []
        b = pysam.AlignedSegment.fromstring(a.to_string(), a.header)
        assert b == a

    def test_string_export_import_str_with_tags(self):
        a = self.build_read()
        a.tags = [("XD", 12), ("RF", "abc")]
        b = pysam.AlignedSegment.fromstring(a.to_string(), a.header)
        assert b == a

    def test_to_string_without_alignment_file(self):
        with open(os.path.join(BAM_DATADIR, "ex2.sam")) as samf:
            reference = [x[:-1] for x in samf if not x.startswith("@")]

        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex2.bam"), "r") as pysamf:
            for s, p in zip(reference, pysamf):
                assert p.to_string() == s

    def test_dict_export(self):
        a = self.build_read()
        a.tags = [("XD", 12), ("RF", "abc")]

        assert a.to_dict() == \
            json.loads(
                '{"name": "read_12345", "flag": "0", "ref_name": "chr1", "ref_pos": "21", '
                '"map_quality": "20", "cigar": "10M1D9M1I20M", "next_ref_name": "=", '
                '"next_ref_pos": "201", "length": "167", '
                '"seq": "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC", '
                '"qual": "1234123412341234123412341234123412341234", "tags": ["XD:i:12", "RF:Z:abc"]}'
            )

    def test_string_export_import_dict_without_tags(self):
        a = self.build_read()
        a.tags = []
        b = pysam.AlignedSegment.from_dict(a.to_dict(), a.header)
        assert b == a

    def test_string_export_import_dict_with_tags(self):
        a = self.build_read()
        a.tags = [("XD", 12), ("RF", "abc")]
        b = pysam.AlignedSegment.from_dict(a.to_dict(), a.header)
        assert b == a


@pytest.mark.parametrize("qual", [
    "",
    "Q",
    pytest.param("""!"#$%&'()*+,-./012...xyz{|}~""", id="linenoise"),
    ">>?AB",
    "ABDDEFGHIJabcdefghij",
    pytest.param("ACAFFGGFFFJDFJHHJIJIHKGGHKHHIJHHHJ7123" * 50, id="long1"),
])
def test_array_to_qualstr(qual):
    qual_array = pysam.qualitystring_to_array(qual)
    result = pysam.array_to_qualitystring(qual_array)
    assert result == qual

def test_longarray_to_qualstr():
    qual_array = array.array('l', [64, 65, 66, 67, 68])
    with pytest.raises(ValueError):
        pysam.array_to_qualitystring(qual_array)


@pytest.mark.parametrize("seq", [
    "A",
    "AT",
    "GCA",
    "ATGC",
    "AATTG",
    pytest.param("ABCDGHKMNRSTVWY", id="iupac"),
])
def test_sequence_unpacking(seq):
    a = pysam.AlignedSegment()
    a.query_sequence = seq
    assert a.query_sequence == seq


@pytest.mark.parametrize("start,stop", [
    (0, 0), (1, 0), (2, 0), (3, 0),
    (0, 1), (1, 1), (2, 1), (3, 1),
    (0, 2), (1, 2), (2, 2), (3, 2),
    (0, 3), (1, 3), (2, 3), (3, 3),
])
def test_subsequence_unpacking(start, stop):
    a = pysam.AlignedSegment()
    a.query_sequence = "ABCDGHKMNRSTVWY"
    start_clip = f"{start}S" if start > 0 else ""
    stop_clip = f"{stop}S" if stop > 0 else ""
    a.cigarstring = f"{start_clip}{a.query_length - start - stop}M{stop_clip}"
    assert a.query_alignment_sequence == a.query_sequence[start : a.query_length - stop]


@pytest.mark.parametrize("seq,revcomp", [
    ("", ""),
    ("A", "T"),
    ("gC", "Gc"),
    ("AAT", "ATT"),
    ("ACGT", "ACGT"),
    ("AATCG", "CGATT"),
    ("aATGGC", "GCCATt"),
    pytest.param("ABCDGHKMNRSTVWY-NNN-abcdghkmnrstvwy--AAASW",
                 "TVGHCDMKNYSABWR-NNN-tvghcdmknysabwr--TTTSW"[::-1], id="iupac"),
])
def test_reverse_complement(seq, revcomp):
    assert isinstance(seq, str)
    assert pysam.reverse_complement(seq) == revcomp

    seq_4byte_kind = f"→{seq}𐘂"
    assert isinstance(seq_4byte_kind, str)
    assert pysam.reverse_complement(seq_4byte_kind) == f"𐘂{revcomp}→"

    seq_bytes = seq.encode("ascii")
    revcomp_bytes = revcomp.encode("ascii")
    assert isinstance(seq_bytes, bytes)
    assert pysam.reverse_complement(seq_bytes) == revcomp_bytes

    seq_ba = bytearray(seq_bytes)
    assert isinstance(seq_ba, bytearray)
    assert pysam.reverse_complement(seq_ba) == revcomp_bytes

    pysam.reverse_complement_inplace(seq_ba)
    assert seq_ba == revcomp_bytes
