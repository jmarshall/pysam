import pysam
import os
import re
import copy
import gzip
from TestUtils import load_and_convert, make_data_files, TABIX_DATADIR


def setUpModule():
    make_data_files(TABIX_DATADIR)


class TestBED:

    filename = os.path.join(TABIX_DATADIR, "fivecolumns.bed.gz")

    def setup_method(self):
        self.tabix = pysam.TabixFile(self.filename)

    def teardown_method(self):
        self.tabix.close()

    def testAssignmentToTargetList(self):
        for row in self.tabix.fetch(parser=pysam.asTuple()):
            # Test that *others gets the right columns...
            contig, start, end, *others = row
            assert 3 + len(others) == len(row)

            # ...and that a TupleProxy can be assigned from more than once
            contig, *others = row
            assert 1 + len(others) == len(row)


class TestParser:

    filename = os.path.join(TABIX_DATADIR, "example.gtf.gz")

    def setup_method(self):
        self.tabix = pysam.TabixFile(self.filename)
        self.compare = load_and_convert(self.filename)

    def teardown_method(self):
        self.tabix.close()

    def testRead(self):
        for x, r in enumerate(self.tabix.fetch(parser=pysam.asTuple())):
            c = self.compare[x]
            assert c == list(r)
            assert len(c) == len(r)

            # test indexing
            for y in range(0, len(r)):
                assert c[y] == r[y]

            # test slicing access
            for y in range(0, len(r) - 1):
                for cc in range(y + 1, len(r)):
                    assert c[y:cc] == r[y:cc]
            assert "\t".join(map(str, c)) == str(r)

    def testAssignmentToTargetList(self):
        for x, r in enumerate(self.tabix.fetch(parser=pysam.asTuple())):
            col1, col2, *others, colN = r
            assert 2 + len(others) + 1 == len(r)

    def testWrite(self):
        for x, r in enumerate(self.tabix.fetch(parser=pysam.asTuple())):
            assert self.compare[x] == list(r)
            c = list(r)
            for y in range(len(r)):
                r[y] = "test_%05i" % y
                c[y] = "test_%05i" % y
            assert [x for x in c] == list(r)
            assert "\t".join(c) == str(r)
            # check second assignment
            for y in range(len(r)):
                r[y] = "test_%05i" % y
            assert [x for x in c] == list(r)
            assert "\t".join(c) == str(r)

    def testUnset(self):
        for x, r in enumerate(self.tabix.fetch(parser=pysam.asTuple())):
            assert self.compare[x] == list(r)
            c = list(r)
            e = list(r)
            for y in range(len(r)):
                r[y] = None
                c[y] = None
                e[y] = ""
                assert c == list(r)
                assert "\t".join(e) == str(r)

    def testIteratorCompressed(self):
        '''test iteration from compressed file.'''
        with gzip.open(self.filename) as infile:
            for x, r in enumerate(pysam.tabix_iterator(
                    infile, pysam.asTuple())):
                assert self.compare[x] == list(r)
                assert len(self.compare[x]) == len(r)

                # test indexing
                for c in range(0, len(r)):
                    assert self.compare[x][c] == r[c]

                # test slicing access
                for c in range(0, len(r) - 1):
                    for cc in range(c + 1, len(r)):
                        assert self.compare[x][c:cc] == r[c:cc]

    def testIteratorUncompressed(self, tmp_path):
        '''test iteration from uncompressed file.'''
        tmpfilename = str(tmp_path / "IteratorUncompressed")
        with gzip.open(self.filename, "rb") as infile, \
                open(tmpfilename, "wb") as outfile:
            outfile.write(infile.read())

        with open(tmpfilename) as infile:
            for x, r in enumerate(pysam.tabix_iterator(
                    infile, pysam.asTuple())):
                assert self.compare[x] == list(r)
                assert len(self.compare[x]) == len(r)

                # test indexing
                for c in range(0, len(r)):
                    assert self.compare[x][c] == r[c]

                # test slicing access
                for c in range(0, len(r) - 1):
                    for cc in range(c + 1, len(r)):
                        assert self.compare[x][c:cc] == r[c:cc]

    def testCopy(self):
        a = next(self.tabix.fetch(parser=pysam.asTuple()))
        b = copy.copy(a)
        assert a == b

        a = next(self.tabix.fetch(parser=pysam.asGTF()))
        b = copy.copy(a)
        assert a == b


class TestGTF(TestParser):

    parser = pysam.asGTF

    def build_attribute_string(self, d):
        """build attribute string from dictionary d"""
        s = "; ".join([f'{x} "{y}"' for (x, y) in d.items()]) + ";"
        # remove quotes around numeric values
        s = re.sub(r'"(\d+)"', r'\1', s)
        return s

    def testRead(self):
        for x, r in enumerate(self.tabix.fetch(parser=self.parser())):
            c = self.compare[x]
            assert len(c) == len(r)
            assert list(c) == list(r)
            assert c == str(r).split("\t")
            assert r.gene_id.startswith("ENSG")
            if r.feature != 'gene':
                assert r.transcript_id.startswith("ENST")
            assert c[0] == r.contig
            assert "\t".join(map(str, c)) == str(r)

    def test_setting_fields(self):
        r = next(self.tabix.fetch(parser=self.parser()))

        r.contig = r.contig + "_test_contig"
        r.source = r.source + "_test_source"
        r.feature = r.feature + "_test_feature"
        r.start += 10
        r.end += 10
        r.score = 20
        r.strand = "+"
        r.frame = 0
        r.attributes = 'gene_id "0001";'
        r.transcript_id = "0002"
        sr = str(r)
        assert "_test_contig" in sr
        assert "_test_source" in sr
        assert "_test_feature" in sr
        assert 'gene_id "0001"' in sr
        assert 'transcript_id "0002"' in sr

    def test_setAttribute_makes_changes(self):
        r = next(self.tabix.fetch(parser=self.parser()))
        r.setAttribute("transcript_id", "abcd")
        sr = str(r)
        assert r.transcript_id == "abcd"
        assert 'transcript_id "abcd"' in sr

    def test_added_attribute_is_output(self):
        r = next(self.tabix.fetch(parser=self.parser()))

        r.new_int_attribute = 12
        assert "new_int_attribute 12" in str(r).split("\t")[8]

        r.new_float_attribute = 12.0
        assert "new_float_attribute 12.0" in str(r).split("\t")[8]

        r.new_text_attribute = "abc"
        assert 'new_text_attribute "abc"' in str(r).split("\t")[8]

    def test_setting_start_is_one_based(self):
        r = next(self.tabix.fetch(parser=self.parser()))
        r.start = 1800
        assert r.start == 1800
        assert str(r).split("\t")[3] == "1801"

    def test_setting_end_is_one_based(self):
        r = next(self.tabix.fetch(parser=self.parser()))
        r.end = 2100
        assert r.end == 2100
        assert str(r).split("\t")[4] == "2100"

    def test_setting_frame_to_none_produces_dot(self):
        r = next(self.tabix.fetch(parser=self.parser()))
        r.frame = None
        assert str(r).split("\t")[7] == "."

        r.frame = 2
        assert str(r).split("\t")[7] == "2"

        r = next(self.tabix.fetch(parser=self.parser()))
        r.frame = "."
        assert r.frame is None
        assert str(r).split("\t")[7] == "."

    def test_setting_source_to_none_produces_dot(self):
        r = next(self.tabix.fetch(parser=self.parser()))
        r.source = None
        assert str(r).split("\t")[1] == "."

        r.source = "source"
        assert str(r).split("\t")[1] == "source"

        r = next(self.tabix.fetch(parser=self.parser()))
        r.source = "."
        assert r.source is None
        assert str(r).split("\t")[1] == "."

    def test_setting_feature_to_none_produces_dot(self):
        r = next(self.tabix.fetch(parser=self.parser()))
        r.feature = None
        assert str(r).split("\t")[2] == "."

        r.feature = "feature"
        assert str(r).split("\t")[2] == "feature"

        r = next(self.tabix.fetch(parser=self.parser()))
        r.feature = "."
        assert r.feature is None
        assert str(r).split("\t")[2] == "."

    def test_setting_strand_to_none_produces_dot(self):
        r = next(self.tabix.fetch(parser=self.parser()))
        r.strand = None
        assert str(r).split("\t")[6] == "."

        r.strand = "-"
        assert str(r).split("\t")[6] == "-"

        r = next(self.tabix.fetch(parser=self.parser()))
        r.strand = "."
        assert r.strand is None
        assert str(r).split("\t")[6] == "."

    def test_setting_score_to_none_produces_dot(self):
        r = next(self.tabix.fetch(parser=self.parser()))
        r.score = None
        assert str(r).split("\t")[5] == "."

        r.score = 12.0
        assert str(r).split("\t")[5] == "12.0"

        r.score = -12.0
        assert str(r).split("\t")[5] == "-12.0"

        r = next(self.tabix.fetch(parser=self.parser()))
        r.score = "."
        assert r.score is None
        assert str(r).split("\t")[5] == "."

        r.score = 12
        assert str(r).split("\t")[5] == "12"

        r.score = -12
        assert str(r).split("\t")[5] == "-12"

    def test_asdict_contains_attributes(self):
        r = next(self.tabix.fetch(parser=self.parser()))
        d = r.to_dict()
        c = self.compare[0]
        s = self.build_attribute_string(d)
        assert s == c[8]

    def test_asdict_can_be_modified(self):
        r = next(self.tabix.fetch(parser=self.parser()))
        d = r.to_dict()
        d["gene_id"] = "new_gene_id"
        expected = 'gene_id "new_gene_id"' if self.parser == pysam.asGTF else "gene_id=new_gene_id"
        assert expected in str(r)


class TestGFF3(TestGTF):

    parser = pysam.asGFF3
    filename = os.path.join(TABIX_DATADIR, "example.gff3.gz")

    def build_attribute_string(self, d):
        """build attribute string from dictionary d"""
        s = ";".join(["{}={}".format(x, y) for (x, y) in d.items()]) + ";"
        return s

    def testRead(self):
        for x, r in enumerate(self.tabix.fetch(parser=self.parser())):
            c = self.compare[x]
            assert len(c) == len(r)
            assert list(c) == list(r)
            assert c == str(r).split("\t")
            assert c[0] == r.contig
            assert "\t".join(map(str, c)) == str(r)
            assert r.ID.startswith("MI00")

    def test_setting_fields(self):
        for r in self.tabix.fetch(parser=self.parser()):
            r.contig = r.contig + "_test_contig"
            r.source = "test_source"
            r.feature = "test_feature"
            r.start += 10
            r.end += 10
            r.score = 20
            r.strand = "+"
            r.frame = 0
            r.ID = "test"
            sr = str(r)
            assert "test_contig" in sr
            assert "test_source" in sr
            assert "test_feature" in sr
            assert "ID=test" in sr

    def test_setAttribute_makes_changes(self):
        r = next(self.tabix.fetch(parser=self.parser()))
        r.setAttribute("transcript_id", "abcd")
        sr = str(r)
        assert r.transcript_id == "abcd"
        assert "transcript_id=abcd" in sr

    def test_added_attribute_is_output(self):
        r = next(self.tabix.fetch(parser=self.parser()))

        r.new_int_attribute = 12
        assert "new_int_attribute=12" in str(r).split("\t")[8]

        r.new_float_attribute = 12.0
        assert "new_float_attribute=12.0" in str(r).split("\t")[8]

        r.new_text_attribute = "abc"
        assert "new_text_attribute=abc" in str(r).split("\t")[8]
