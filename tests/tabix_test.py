import os
import shutil
import gzip
import pysam
import glob
import pytest
import re

from TestUtils import checkBinaryEqual, checkGZBinaryEqual, check_url, \
    load_and_convert, make_data_files, TABIX_DATADIR


def setUpModule():
    make_data_files(TABIX_DATADIR)


def myzip_open(infile, mode="r"):
    '''open compressed file and decode.'''

    def _convert(f):
        for l in f:
            yield l.decode("ascii")

    if mode == "r":
        return _convert(gzip.open(infile, "r"))


def splitToBytes(s):
    '''split string and return list of bytes.'''
    return [x.encode("ascii") for x in s.split("\t")]


class TestIndexing:
    filename = os.path.join(TABIX_DATADIR, "example.gtf.gz")
    filename_idx = os.path.join(TABIX_DATADIR, "example.gtf.gz.tbi")

    @pytest.fixture(autouse=True)
    def copy(self, tmp_path):
        self.tmpfilename = str(tmp_path / "copy.gtf.gz")
        shutil.copyfile(self.filename, self.tmpfilename)

    def test_indexing_with_preset_works(self):
        '''test indexing via preset.'''

        pysam.tabix_index(self.tmpfilename, preset="gff")
        assert checkGZBinaryEqual(self.tmpfilename + ".tbi", self.filename_idx)

    def test_indexing_to_custom_location_works(self, tmp_path):
        '''test indexing a file with a non-default location.'''

        index_path = str(tmp_path / "custom.tbi")
        pysam.tabix_index(self.tmpfilename, preset="gff",
                          index=index_path, force=True)
        assert checkGZBinaryEqual(index_path, self.filename_idx)

    def test_indexing_with_explict_columns_works(self):
        '''test indexing via preset.'''

        pysam.tabix_index(self.tmpfilename,
                          seq_col=0,
                          start_col=3,
                          end_col=4,
                          line_skip=0,
                          zerobased=False)
        assert checkGZBinaryEqual(self.tmpfilename + ".tbi", self.filename_idx)

    def test_indexing_with_lineskipping_works(self):
        '''test indexing via preset and lineskip.'''
        pysam.tabix_index(self.tmpfilename,
                          seq_col=0,
                          start_col=3,
                          end_col=4,
                          line_skip=1,
                          zerobased=False)
        assert not checkGZBinaryEqual(self.tmpfilename + ".tbi", self.filename_idx)


class TestCompression:
    filename = os.path.join(TABIX_DATADIR, "example.gtf.gz")
    filename_idx = os.path.join(TABIX_DATADIR, "example.gtf.gz.tbi")
    preset = "gff"

    @pytest.fixture(autouse=True)
    def copy(self, tmp_path):
        self.tmpfilename = str(tmp_path / "copy.gtf")
        with gzip.open(self.filename, "rb") as infile, \
                open(self.tmpfilename, "wb") as outfile:
            outfile.write(infile.read())

    def testCompression(self):
        '''see also issue 106'''
        pysam.tabix_compress(self.tmpfilename, self.tmpfilename + ".gz")
        checkBinaryEqual(self.tmpfilename, self.tmpfilename + ".gz")

    def testIndexPresetUncompressed(self):
        '''test indexing via preset.'''

        pysam.tabix_index(self.tmpfilename, preset=self.preset)
        # check if uncompressed file has been removed
        assert not os.path.exists(self.tmpfilename)
        checkBinaryEqual(self.tmpfilename + ".gz", self.filename)
        checkBinaryEqual(self.tmpfilename + ".gz.tbi", self.filename_idx)

    def testIndexPresetCompressed(self):
        '''test indexing via preset.'''

        pysam.tabix_compress(self.tmpfilename, self.tmpfilename + ".gz")
        pysam.tabix_index(self.tmpfilename + ".gz", preset=self.preset)
        checkBinaryEqual(self.tmpfilename + ".gz", self.filename)
        checkBinaryEqual(self.tmpfilename + ".gz.tbi", self.filename_idx)


class TestCompressionSam(TestCompression):
    filename = os.path.join(TABIX_DATADIR, "example.sam.gz")
    filename_index = os.path.join(TABIX_DATADIR, "example.sam.gz.tbi")
    preset = "sam"


class TestCompressionBed(TestCompression):
    filename = os.path.join(TABIX_DATADIR, "example.bed.gz")
    filename_index = os.path.join(TABIX_DATADIR, "example.bed.gz.tbi")
    preset = "bed"


class TestCompressionVCF(TestCompression):
    filename = os.path.join(TABIX_DATADIR, "example.vcf.gz")
    filename_index = os.path.join(TABIX_DATADIR, "example.vcf.gz.tbi")
    preset = "vcf"


class IterationTest:

    with_comments = False

    def setup_method(self):
        lines = []
        with gzip.open(self.filename, "rb") as inf:
            for line in inf:
                line = line.decode('ascii')
                if line.startswith("#"):
                    if not self.with_comments:
                        continue
                lines.append(line)

        # creates index of contig, start, end, adds content without newline.
        self.compare = [
            (x[0][0], int(x[0][3]), int(x[0][4]), x[1])
            for x in [(y.split("\t"), y[:-1]) for y in lines
                      if not y.startswith("#")]]

        self.comments = [x[:-1] for x in lines if x.startswith("#")]

    def getSubset(self, contig=None, start=None, end=None):

        if contig is None:
            # all lines
            subset = [x[3] for x in self.compare]
        else:
            if start is not None and end is None:
                # until end of contig
                subset = [x[3]
                          for x in self.compare if x[0] == contig and
                          x[2] > start]
            elif start is None and end is not None:
                # from start of contig
                subset = [x[3]
                          for x in self.compare if x[0] == contig and
                          x[1] <= end]
            elif start is None and end is None:
                subset = [x[3] for x in self.compare if x[0] == contig]
            else:
                # all within interval
                subset = [x[3] for x in self.compare if x[0] == contig and
                          min(x[2], end) - max(x[1], start) > 0]

        if self.with_comments:
            subset.extend(self.comments)

        return subset

    def checkPairwise(self, result, ref):
        '''check pairwise results.
        '''
        result.sort()
        ref.sort()
        a = set(result)
        b = set(ref)

        assert len(result) == len(ref), f"differences are {a.difference(b)}, {b.difference(a)}"

        for x, d in enumerate(list(zip(result, ref))):
            assert d[0] == d[1], f"unexpected results in pair {x}:\n{d[0]!r}', expected\n{d[1]!r}"


class TestGZFile(IterationTest):

    filename = os.path.join(TABIX_DATADIR, "example.gtf.gz")
    with_comments = True

    def setup_method(self):
        super().setup_method()
        self.gzfile = pysam.GZIterator(self.filename)

    def testAll(self):
        result = list(self.gzfile)
        ref = self.getSubset()
        self.checkPairwise(result, ref)


class TestIterationWithoutComments(IterationTest):

    '''test iterating with TabixFile.fetch() when
    there are no comments in the file.'''

    filename = os.path.join(TABIX_DATADIR,
                            "example.gtf.gz")

    def setup_method(self):
        super().setup_method()
        self.tabix = pysam.TabixFile(self.filename)

    def teardown_method(self):
        self.tabix.close()

    def testRegionStrings(self):
        """test if access with various region strings works"""

        assert len(list(self.tabix.fetch("chr1"))) == 218
        assert len(list(self.tabix.fetch("chr1", 1000))) == 218
        assert len(list(self.tabix.fetch("chr1", end=1000000))) == 218
        assert len(list(self.tabix.fetch("chr1", 1000, 1000000))) == 218

    def testAll(self):
        result = list(self.tabix.fetch())
        ref = self.getSubset()
        self.checkPairwise(result, ref)

    def testPerContig(self):
        for contig in ("chr1", "chr2", "chr1", "chr2"):
            result = list(self.tabix.fetch(contig))
            ref = self.getSubset(contig)
            self.checkPairwise(result, ref)

    def testPerContigToEnd(self):

        end = None
        for contig in ("chr1", "chr2", "chr1", "chr2"):
            for start in range(0, 200000, 1000):
                result = list(self.tabix.fetch(contig, start, end))
                ref = self.getSubset(contig, start, end)
                self.checkPairwise(result, ref)

    def testPerContigFromStart(self):

        start = None
        for contig in ("chr1", "chr2", "chr1", "chr2"):
            for end in range(0, 200000, 1000):
                result = list(self.tabix.fetch(contig, start, end))
                ref = self.getSubset(contig, start, end)
                self.checkPairwise(result, ref)

    def testPerContig2(self):

        start, end = None, None
        for contig in ("chr1", "chr2", "chr1", "chr2"):
            result = list(self.tabix.fetch(contig, start, end))
            ref = self.getSubset(contig, start, end)
            self.checkPairwise(result, ref)

    def testPerInterval(self):

        start, end = None, None
        for contig in ("chr1", "chr2", "chr1", "chr2"):
            for start in range(0, 200000, 2000):
                for end in range(start, start + 2000, 500):
                    result = list(self.tabix.fetch(contig, start, end))
                    ref = self.getSubset(contig, start, end)
                    self.checkPairwise(result, ref)

    def testInvalidIntervals(self):

        # invalid intervals (start > end)
        with pytest.raises(ValueError): self.tabix.fetch("chr1", 0, -10)
        with pytest.raises(ValueError): self.tabix.fetch("chr1", 200, 0)

        # out of range intervals
        with pytest.raises(ValueError): self.tabix.fetch("chr1", -10, 200)
        with pytest.raises(ValueError): self.tabix.fetch("chr1", -10, -20)

        # unknown chromosome
        with pytest.raises(ValueError): self.tabix.fetch("chrUn")

        # out of range access
        # to be implemented
        # wth pytest.raises(IndexError): self.tabix.fetch("chr1", 1000000, 2000000)

        # raise no error for empty intervals
        self.tabix.fetch("chr1", 100, 100)

    def testGetContigs(self):
        assert sorted(self.tabix.contigs) == ["chr1", "chr2"]
        # check that contigs is read-only
        with pytest.raises(AttributeError): self.tabix.contigs = ["chr1", "chr2"]

    def testHeader(self):
        ref = []
        with gzip.open(self.filename) as inf:
            for x in inf:
                x = x.decode("ascii")
                if not x.startswith("#"):
                    break
                ref.append(x[:-1])

        header = list(self.tabix.header)
        assert ref == header

    def testReopening(self):
        '''test repeated opening of the same file.'''
        def func1():
            # opens any tabix file
            with pysam.TabixFile(self.filename) as inf:
                pass

        for i in range(1000):
            func1()


class TestIterationWithComments(TestIterationWithoutComments):

    '''test iterating with TabixFile.fetch() when
    there are comments in the file.

    Tests will create plenty of warnings on stderr.
    '''

    filename = os.path.join(TABIX_DATADIR, "example_comments.gtf.gz")


class TestIterators:
    filename = os.path.join(TABIX_DATADIR, "example.gtf.gz")

    iterator = pysam.tabix_generic_iterator
    parser = pysam.asTuple
    is_compressed = False

    @pytest.fixture(autouse=True)
    def copy(self, tmp_path):
        self.tabix = pysam.TabixFile(self.filename)
        self.compare = load_and_convert(self.filename)
        self.tmpfilename_uncompressed = str(tmp_path / "TestIterators")
        with gzip.open(self.filename, "rb") as infile, \
                open(self.tmpfilename_uncompressed, "wb") as outfile:
            outfile.write(infile.read())

        yield

        self.tabix.close()

    def open(self):
        if self.is_compressed:
            infile = gzip.open(self.filename)
        else:
            infile = open(self.tmpfilename_uncompressed)
        return infile

    def testIteration(self):
        with self.open() as infile:
            for x, r in enumerate(self.iterator(infile, self.parser())):
                assert self.compare[x] == list(r)
                assert len(self.compare[x]) == len(r)

                # test indexing
                for c in range(0, len(r)):
                    assert self.compare[x][c] == r[c]

                # test slicing access
                for c in range(0, len(r) - 1):
                    for cc in range(c + 1, len(r)):
                        assert self.compare[x][c:cc] == r[c:cc]

    def testClosedFile(self):
        '''test for error when iterating from closed file.'''
        infile = self.open()
        infile.close()

        # iterating from a closed file should raise a value error
        with pytest.raises(ValueError): self.iterator(infile, self.parser())

    def testClosedFileIteration(self):
        '''test for error when iterating from file that has been closed'''

        infile = self.open()

        i = self.iterator(infile, self.parser())
        x = next(i)
        infile.close()
        # Not implemented
        # with pytest.raises(ValueError): next(i)


class TestIteratorsGenericCompressed(TestIterators):
    is_compressed = True


class TestIteratorsFileCompressed(TestIterators):
    iterator = pysam.tabix_file_iterator
    is_compressed = True


class TestIteratorsFileUncompressed(TestIterators):
    iterator = pysam.tabix_file_iterator
    is_compressed = False


class TestIterationMalformattedGTFFiles:

    '''test reading from malformatted gtf files.'''

    iterator = pysam.tabix_generic_iterator
    parser = pysam.asGTF

    def testGTFTooManyFields(self):

        with gzip.open(os.path.join(
                TABIX_DATADIR,
                "gtf_toomany_fields.gtf.gz")) as infile:
            iterator = self.iterator(
                infile,
                parser=self.parser())
            with pytest.raises(ValueError): next(iterator)

    def testGTFTooFewFields(self):

        with gzip.open(os.path.join(
                TABIX_DATADIR,
                "gtf_toofew_fields.gtf.gz")) as infile:
            iterator = self.iterator(
                infile,
                parser=self.parser())
            with pytest.raises(ValueError): next(iterator)


class TestBed:
    filename = os.path.join(TABIX_DATADIR, "example.bed.gz")

    def setup_method(self):
        self.tabix = pysam.TabixFile(self.filename)
        self.compare = load_and_convert(self.filename)

    def teardown_method(self):
        self.tabix.close()

    def testRead(self):

        for x, r in enumerate(self.tabix.fetch(parser=pysam.asBed())):
            c = self.compare[x]
            assert len(c) == len(r)
            assert c == str(r).split("\t")
            assert c[0] == r.contig
            assert int(c[1]) == r.start
            assert int(c[2]) == r.end
            # Needs lambda so that the property getter isn't called too soon
            with pytest.raises(KeyError): r.name
            with pytest.raises(KeyError): r.score
            assert list(c) == list(r)
            assert "\t".join(map(str, c)) == str(r)

    def testWrite(self):

        for x, r in enumerate(self.tabix.fetch(parser=pysam.asBed())):
            c = self.compare[x]
            assert c == str(r).split("\t")
            assert list(c) == list(r)

            r.contig = "test"
            assert r.contig == "test"
            assert r[0] == "test"

            r.start += 1
            assert int(c[1]) + 1 == r.start
            assert str(int(c[1]) + 1) == r[1]

            r.end += 1
            assert int(c[2]) + 1 == r.end
            assert str(int(c[2]) + 1) == r[2]

            with pytest.raises(IndexError):
                r.name = "test"

            with pytest.raises(IndexError):
                r.score = 1


class TestVCF:

    filename = os.path.join(TABIX_DATADIR, "example.vcf40")

    @pytest.fixture(autouse=True)
    def copy_and_index(self, tmp_path):
        self.tmpfilename = str(tmp_path / "copy.vcf")
        shutil.copyfile(self.filename, self.tmpfilename)
        pysam.tabix_index(self.tmpfilename, preset="vcf")


class TestUnicode:

    '''test reading from a file with non-ascii characters.'''

    filename = os.path.join(TABIX_DATADIR, "example_unicode.vcf")

    @pytest.fixture(autouse=True)
    def copy_and_index(self, tmp_path):
        self.tmpfilename = str(tmp_path / "copy.vcf")
        shutil.copyfile(self.filename, self.tmpfilename)
        pysam.tabix_index(self.tmpfilename, preset="vcf")

    def testFromTabix(self):
        # use ascii encoding - should raise error
        with pysam.TabixFile(
                self.tmpfilename + ".gz", encoding="ascii") as t:
            results = list(t.fetch(parser=pysam.asVCF()))
            with pytest.raises(UnicodeDecodeError):
                results[1].id

        with pysam.TabixFile(
                self.tmpfilename + ".gz", encoding="utf-8") as t:
            results = list(t.fetch(parser=pysam.asVCF()))
            assert results[1].id == "Reneé"

    def testFromVCF(self):
        self.vcf = pysam.VCF()
        with pytest.raises(UnicodeDecodeError):
            self.vcf.connect(self.tmpfilename + ".gz", "ascii")
        self.vcf.connect(self.tmpfilename + ".gz", encoding="utf-8")
        v = self.vcf.getsamples()[0]


class TestVCFFromTabix(TestVCF):

    columns = ("contig", "pos", "id",
               "ref", "alt", "qual",
               "filter", "info", "format")

    @pytest.fixture(autouse=True)
    def tabix_and_load(self):
        self.tabix = pysam.TabixFile(self.tmpfilename + ".gz")
        self.compare = load_and_convert(self.filename)

        yield

        self.tabix.close()

    def testRead(self):
        ncolumns = len(self.columns)

        for x, r in enumerate(self.tabix.fetch(parser=pysam.asVCF())):
            c = self.compare[x]
            for y, field in enumerate(self.columns):
                # it is ok to have a missing format column
                if y == 8 and y == len(c):
                    continue
                if field == "pos":
                    assert int(c[y]) - 1 == getattr(r, field)
                    assert int(c[y]) - 1 == r.pos
                else:
                    assert c[y] == getattr(r, field), f"mismatch in {field}"
            if len(c) == 8:
                assert len(r) == 0
            else:
                assert len(c) == len(r) + ncolumns

            for y in range(len(c) - ncolumns):
                assert c[ncolumns + y] == r[y]
            assert "\t".join(map(str, c)) == str(r)

    def testWrite(self):
        ncolumns = len(self.columns)

        for x, r in enumerate(self.tabix.fetch(parser=pysam.asVCF())):
            c = self.compare[x]
            # check unmodified string
            cmp_string = str(r)
            ref_string = "\t".join([x for x in c])

            assert ref_string == cmp_string

            # set fields and compare field-wise
            for y, field in enumerate(self.columns):
                # it is ok to have a missing format column
                if y == 8 and y == len(c):
                    continue
                if field == "pos":
                    rpos = getattr(r, field)
                    assert int(c[y]) - 1 == rpos
                    assert int(c[y]) - 1 == r.pos
                    # increment pos by 1
                    setattr(r, field, rpos + 1)
                    assert getattr(r, field) == rpos + 1
                    c[y] = str(int(c[y]) + 1)
                else:
                    setattr(r, field, "test_%i" % y)
                    c[y] = "test_%i" % y
                    assert c[y] == getattr(r, field), f"mismatch in field {field}"

            if len(c) == 8:
                assert len(r) == 0
            else:
                assert len(c) == len(r) + ncolumns

            for y in range(len(c) - ncolumns):
                c[ncolumns + y] = "test_%i" % y
                r[y] = "test_%i" % y
                assert c[ncolumns + y] == r[y]


class TestVCFFromVCF(TestVCF):

    columns = ("chrom", "pos", "id",
               "ref", "alt", "qual",
               "filter", "info", "format")

    # tests failing while parsing
    fail_on_parsing = (
        (5, "Flag fields should not have a value"),
        (9, "aouao"),
        (12, "Error BAD_NUMBER_OF_PARAMETERS"),
        (13, "aoeu"),
        (18, "Error BAD_NUMBER_OF_PARAMETERS"),
        (24, "Error HEADING_NOT_SEPARATED_BY_TABS"))

    # tests failing on opening
    fail_on_opening = ((24, "Error HEADING_NOT_SEPARATED_BY_TABS"),
                       )

    fail_on_samples = []

    check_samples = False
    coordinate_offset = 1

    # value returned for missing values
    missing_value = "."
    missing_quality = -1

    @pytest.fixture(autouse=True)
    def vcf_and_load(self):
        self.vcf = pysam.VCF()
        self.compare = load_and_convert(self.filename, encode=False)

        yield

        self.vcf.close()

    def open_vcf(self, fn):
        return self.vcf.connect(fn)

    def get_failure_stage(self):
        fn = os.path.basename(self.filename)
        for x, msg in self.fail_on_opening:
            if "{}.vcf".format(x) == fn:
                return "opening"

        for x, msg in self.fail_on_parsing:
            if "{}.vcf".format(x) == fn:
                return "parsing"

        for x, msg in self.fail_on_samples:
            if "{}.vcf".format(x) == fn:
                return "samples"

        return None

    def testConnecting(self):
        if self.get_failure_stage() == "opening":
            with pytest.raises(ValueError):
                self.open_vcf(self.tmpfilename + ".gz")
        else:
            self.open_vcf(self.tmpfilename + ".gz")

    def get_iterator(self):
        with open(self.filename) as f:
            fn = os.path.basename(self.filename)
            return list(self.vcf.parse(f))

    def get_field_value(self, record, field):
        return record[field]

    def sample2value(self, r, v):
        return r, v

    def alt2value(self, r, v):
        if r == ".":
            return [], v
        else:
            return r.split(","), list(v)

    def filter2value(self, r, v):
        if r == "PASS":
            return [], v
        elif r == ".":
            return [], v
        else:
            return r.split(";"), v

    def testParsing(self):
        if self.get_failure_stage() in ("opening", "parsing"):
            return

        itr = self.get_iterator()
        if itr is None:
            return

        fn = os.path.basename(self.filename)

        check_samples = self.check_samples
        for vcf_code, msg in self.fail_on_samples:
            if "%i.vcf" % vcf_code == fn:
                check_samples = False

        for x, r in enumerate(itr):
            c = self.compare[x]
            for y, field in enumerate(self.columns):
                # it is ok to have a missing format column
                if y == 8 and y == len(c):
                    continue

                val = self.get_field_value(r, field)
                if field == "pos":
                    assert int(c[y]) - self.coordinate_offset == val
                elif field == "alt" or field == "alts":
                    cc, vv = self.alt2value(c[y], val)
                    if cc != vv:
                        # import pdb; pdb.set_trace()
                        pass
                    assert cc == vv, f"mismatch in {field}"

                elif field == "filter":
                    cc, vv = self.filter2value(c[y], val)
                    assert cc == vv, f"mismatch in {field}"

                elif field == "info":
                    # tests for info field not implemented
                    pass

                elif field == "qual" and c[y] == ".":
                    assert val == self.missing_quality, f"mismatch in {field}"

                elif field == "format":
                    # format field converted to list
                    assert c[y].split(":") == list(val), f"mismatch in {field}"

                elif type(val) in (int, float):
                    if c[y] == ".":
                        assert val is None, f"mismatch in {field}"
                    else:
                        assert float(c[y]) == pytest.approx(float(val)), f"mismatch in {field}"
                else:
                    if c[y] == ".":
                        ref_val = self.missing_value
                    else:
                        ref_val = c[y]
                    assert val == ref_val, f"mismatch in {field}"
            # parse samples
            if check_samples:
                if len(c) == 8:
                    for x, s in enumerate(r.samples):
                        assert r.samples[s].values() == [], f"mismatch in sample {s}"
                else:
                    for x, s in enumerate(r.samples):
                        ref, comp = self.sample2value(
                            c[9 + x],
                            r.samples[s])
                        self.compare_samples(ref, comp, s, r)

    def compare_samples(self, ref, comp, s, r):

        if ref != comp:

            # check if GT not at start, not VCF conform and
            # not supported by cbcf.pyx
            k = r.format.keys()
            if "GT" in k and k[0] != "GT":
                return

            # perform an element-wise checto work around rounding differences
            for a, b in zip(re.split("[:,;]", ref),
                            re.split("[:,;]", comp)):
                is_float = True
                try:
                    a = float(a)
                    b = float(b)
                except ValueError:
                    is_float = False

                if is_float:
                    assert a == pytest.approx(b), "mismatch in sample {s}: expected {ref}, got {comp}"
                else:
                    assert a == b, "mismatch in sample {s}: expected {ref}, got {comp}"


############################################################################
# create a test class for each example vcf file.
# Two samples are created -
# 1. Testing pysam/tabix access
# 2. Testing the VCF class
vcf_files = glob.glob(os.path.join(TABIX_DATADIR, "vcf", "*.vcf"))

for vcf_file in vcf_files:
    n = "TestVCFFromTabix_%s" % os.path.basename(vcf_file[:-4])
    globals()[n] = type(n, (TestVCFFromTabix,), dict(filename=vcf_file,))
    n = "TestVCFFromVCF_%s" % os.path.basename(vcf_file[:-4])
    globals()[n] = type(n, (TestVCFFromVCF,), dict(filename=vcf_file,))


class TestVCFFromVariantFile(TestVCFFromVCF):

    columns = ("chrom", "pos", "id",
               "ref", "alts", "qual",
               "filter", "info", "format")

    fail_on_parsing = [
        (24, 'Could not parse the "#CHROM.." line'),
        ("issue85", "empty VCF"),
    ]
    fail_on_opening = [
        (24, 'Could not parse the "#CHROM.." line'),
        ("issue85", "empty VCF"),
    ]
    coordinate_offset = 0
    check_samples = True
    fail_on_samples = [
        (9, "PL field not defined. Expected to be scalar, but is array"),
        (12, "PL field not defined. Expected to be scalar, but is array"),
        (18, "PL field not defined. Expected to be scalar, but is array"),
    ]

    # value returned for missing values
    missing_value = None
    missing_quality = None

    vcf = None

    def filter2value(self, r, v):
        if r == "PASS":
            return ["PASS"], list(v)
        elif r == ".":
            return [], list(v)
        else:
            return r.split(";"), list(v)

    def alt2value(self, r, v):
        if r == ".":
            return None, v
        else:
            return r.split(","), list(v)

    def sample2value(self, r, smp):

        def convert_field(f):
            if f is None:
                return "."
            elif isinstance(f, tuple):
                return ",".join(map(convert_field, f))
            else:
                return str(f)

        v = smp.values()

        if 'GT' in smp:
            alleles = [
                str(a) if a is not None else '.' for a in smp.allele_indices]
            v[0] = '/|'[smp.phased].join(alleles)

        comp = ":".join(map(convert_field, v))

        if comp.endswith(":."):
            comp = comp[:-2]

        return r, comp

    def setup_method(self):
        self.compare = load_and_convert(self.filename, encode=False)

    def teardown_method(self):
        if self.vcf:
            self.vcf.close()
        self.vcf = None

    def get_iterator(self):
        self.vcf = pysam.VariantFile(self.filename)
        return self.vcf.fetch()

    def get_field_value(self, record, field):
        return getattr(record, field)

    def open_vcf(self, fn):
        with pysam.VariantFile(fn) as inf:
            pass


for vcf_file in vcf_files:
    p = os.path.basename(vcf_file[:-4])
    n = "TestVCFFromVariantFile_%s" % p
    globals()[n] = type(n, (TestVCFFromVariantFile,), dict(filename=vcf_file,))


class TestRemoteFileHTTP:

    url = "http://genserv.anat.ox.ac.uk/downloads/pysam/test/example.gtf.gz"
    region = "chr1:1-1000"
    local = os.path.join(TABIX_DATADIR, "example.gtf.gz")

    def setup_method(self):
        if not getattr(pysam.config, "HAVE_LIBCURL", 0) or not check_url(self.url):
            self.remote_file = None
        else:
            self.remote_file = pysam.TabixFile(self.url, "r")

        self.local_file = pysam.TabixFile(self.local, "r")

    def teardown_method(self):
        if self.remote_file is None:
            return

        self.remote_file.close()
        self.local_file.close()

    def testFetchAll(self):
        if self.remote_file is None:
            return

        remote_result = list(self.remote_file.fetch())
        local_result = list(self.local_file.fetch())

        assert len(remote_result) == len(local_result)
        for x, y in zip(remote_result, local_result):
            assert x == y

    def testHeader(self):
        if self.remote_file is None:
            return

        assert list(self.local_file.header) == []


class TestRemoteFileHTTPWithHeader(TestRemoteFileHTTP):

    url = "http://genserv.anat.ox.ac.uk/downloads/pysam/test/example_comments.gtf.gz"
    region = "chr1:1-1000"
    local = os.path.join(TABIX_DATADIR, "example_comments.gtf.gz")

    def setup_method(self):
        if not getattr(pysam.config, "HAVE_LIBCURL", 0) or not check_url(self.url):
            self.remote_file = None
        else:
            self.remote_file = pysam.TabixFile(self.url, "r")
        self.local_file = pysam.TabixFile(self.local, "r")

    def testHeader(self):
        if self.remote_file is None:
            return

        assert list(self.local_file.header) == ["# comment at start"]
        assert list(self.local_file.header) == self.remote_file.header


class TestIndexArgument:

    filename_src = os.path.join(TABIX_DATADIR, "example.vcf.gz")
    index_src = os.path.join(TABIX_DATADIR, "example.vcf.gz.tbi")

    def testFetchAll(self, tmp_path):
        filename_dst = str(tmp_path / "example.vcf.gz")
        index_dst    = str(tmp_path / "example.vcf.gz.tbi")

        shutil.copyfile(self.filename_src, filename_dst)
        shutil.copyfile(self.index_src, index_dst)

        with pysam.TabixFile(
                self.filename_src, "r", index=self.index_src) as same_basename_file:
            same_basename_results = list(same_basename_file.fetch())

        with pysam.TabixFile(
                filename_dst, "r", index=index_dst) as diff_index_file:
            diff_index_result = list(diff_index_file.fetch())

        assert len(same_basename_results) == len(diff_index_result)
        for x, y in zip(same_basename_results, diff_index_result):
            assert x == y

    def testLoadIndexWithoutTbiExtension(self, tmp_path):
        filename_dst  = str(tmp_path / "example.vcf.gz")
        index_dst_dat = str(tmp_path / "example.vcf.gz.tbi.dat")

        shutil.copyfile(self.filename_src, filename_dst)
        shutil.copyfile(self.index_src, index_dst_dat)

        with pysam.TabixFile(self.filename_src, "r", index=self.index_src) as same_basename_file:
            same_basename_results = list(same_basename_file.fetch())

        with pysam.TabixFile(filename_dst, "r", index=index_dst_dat) as diff_index_file:
            diff_index_result = list(diff_index_file.fetch())

        assert len(same_basename_results) == len(diff_index_result)
        for x, y in zip(same_basename_results, diff_index_result):
            assert x == y


def _TestMultipleIteratorsHelper(filename, multiple_iterators):
    '''open file within scope, return iterator.'''

    tabix = pysam.TabixFile(filename)
    iterator = tabix.fetch(parser=pysam.asGTF(),
                           multiple_iterators=multiple_iterators)
    tabix.close()
    return iterator


class TestBackwardsCompatibility:
    """check if error is raised if a tabix file from an
    old version is accessed from pysam"""

    def check(self, filename, raises=None):
        with pysam.TabixFile(filename) as tf:
            ref = load_and_convert(filename)
            if raises is None:
                assert len(list(tf.fetch())) == len(ref)
            else:
                with pytest.raises(raises): tf.fetch()

    def testVCF0v23(self):
        self.check(os.path.join(TABIX_DATADIR, "example_0v23.vcf.gz"),
                   ValueError)

    def testBED0v23(self):
        self.check(os.path.join(TABIX_DATADIR, "example_0v23.bed.gz"),
                   ValueError)

    def testVCF0v26(self):
        self.check(os.path.join(TABIX_DATADIR, "example_0v26.vcf.gz"),
                   ValueError)

    def testBED0v26(self):
        self.check(os.path.join(TABIX_DATADIR, "example_0v26.bed.gz"),
                   ValueError)

    def testVCF(self):
        self.check(os.path.join(TABIX_DATADIR, "example.vcf.gz"))

    def testBED(self):
        self.check(os.path.join(TABIX_DATADIR, "example.bed.gz"))

    def testEmpty(self):
        self.check(os.path.join(TABIX_DATADIR, "empty.bed.gz"))


class TestMultipleIterators:

    filename = os.path.join(TABIX_DATADIR, "example.gtf.gz")

    def testJoinedIterators(self):

        # two iterators working on the same file
        with pysam.TabixFile(self.filename) as tabix:
            a = next(tabix.fetch(parser=pysam.asGTF()))
            b = next(tabix.fetch(parser=pysam.asGTF()))
            # the first two lines differ only by the feature field
            assert a.feature == "UTR"
            assert b.feature == "exon"
            assert re.sub("UTR", "", str(a)) == re.sub("exon", "", str(b))

    def testDisjointIterators(self):
        # two iterators working on the same file
        with pysam.TabixFile(self.filename) as tabix:
            a = next(tabix.fetch(parser=pysam.asGTF(), multiple_iterators=True))
            b = next(tabix.fetch(parser=pysam.asGTF(), multiple_iterators=True))
            # both iterators are at top of file
            assert str(a) == str(b)

    def testScope(self):
        # technically it does not really test if the scope is correct
        i = _TestMultipleIteratorsHelper(self.filename,
                                         multiple_iterators=True)
        assert next(i)
        i = _TestMultipleIteratorsHelper(self.filename,
                                         multiple_iterators=False)
        with pytest.raises(IOError): next(i)

    def testDoubleFetch(self):
        with pysam.TabixFile(self.filename) as f:
            for a, b in zip(f.fetch(multiple_iterators=True),
                            f.fetch(multiple_iterators=True)):
                assert str(a) == str(b)


class TestContextManager:

    filename = os.path.join(TABIX_DATADIR, "example.gtf.gz")

    def testManager(self):
        with pysam.TabixFile(self.filename) as tabixfile:
            tabixfile.fetch()
        assert tabixfile.closed


class TestMultithreadTabixFile:

    filename = os.path.join(TABIX_DATADIR, "example.gtf.gz")

    def testMultithreadEqualsSinglethread(self):
        with pysam.TabixFile(self.filename) as tabixfile:
            single = [r for r in tabixfile.fetch()]
        with pysam.TabixFile(self.filename, threads=2) as tabixfile:
            multi = [r for r in tabixfile.fetch()]
        for r1, r2 in zip(single, multi):
            assert str(r1) == str(r2)
