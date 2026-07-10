import pysam
import os
import gzip
import copy
import shutil

import pytest

from TestUtils import check_url, make_data_files, BAM_DATADIR, get_temp_filename


def setUpModule():
    make_data_files(BAM_DATADIR)


class TestFastaFile:

    sequences = {
        'chr1':
        "CACTAGTGGCTCATTGTAAATGTGTGGTTTAACTCGTCCATGGCCCAGCATTAGGGAGCTGTGGACCCTGCAGCCTGGCTGTGGGGGCCGCAGTGGCTGAGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCACAGGGGCTTAACCTCTGGTGACTGCCAGAGCTGCTGGCAAGCTAGAGTCCCATTTGGAGCCCCTCTAAGCCGTTCTATTTGTAATGAAAACTATATTTATGCTATTCAGTTCTAAATATAGAAATTGAAACAGCTGTGTTTAGTGCCTTTGTTCAACCCCCTTGCAACAACCTTGAGAACCCCAGGGAATTTGTCAATGTCAGGGAAGGAGCATTTTGTCAGTTACCAAATGTGTTTATTACCAGAGGGATGGAGGGAAGAGGGACGCTGAAGAACTTTGATGCCCTCTTCTTCCAAAGATGAAACGCGTAACTGCGCTCTCATTCACTCCAGCTCCCTGTCACCCAATGGACCTGTGATATCTGGATTCTGGGAAATTCTTCATCCTGGACCCTGAGAGATTCTGCAGCCCAGCTCCAGATTGCTTGTGGTCTGACAGGCTGCAACTGTGAGCCATCACAATGAACAACAGGAAGAAAAGGTCTTTCAAAAGGTGATGTGTGTTCTCATCAACCTCATACACACACATGGTTTAGGGGTATAATACCTCTACATGGCTGATTATGAAAACAATGTTCCCCAGATACCATCCCTGTCTTACTTCCAGCTCCCCAGAGGGAAAGCTTTCAACGCTTCTAGCCATTTCTTTTGGCATTTGCCTTCAGACCCTACACGAATGCGTCTCTACCACAGGGGGCTGCGCGGTTTCCCATCATGAAGCACTGAACTTCCACGTCTCATCTAGGGGAACAGGGAGGTGCACTAATGCGCTCCACGCCCAAGCCCTTCTCACAGTTTCTGCCCCCAGCATGGTTGTACTGGGCAATACATGAGATTATTAGGAAATGCTTTACTGTCATAACTATGAAGAGACTATTGCCAGATGAACCACACATTAATACTATGTTTCTTATCTGCACATTACTACCCTGCAATTAATATAATTGTGTCCATGTACACACGCTGTCCTATGTACTTATCATGACTCTATCCCAAATTCCCAATTACGTCCTATCTTCTTCTTAGGGAAGAACAGCTTAGGTATCAATTTGGTGTTCTGTGTAAAGTCTCAGGGAGCCGTCCGTGTCCTCCCATCTGGCCTCGTCCACACTGGTTCTCTTGAAAGCTTGGGCTGTAATGATGCCCCTTGGCCATCACCCAGTCCCTGCCCCATCTCTTGTAATCTCTCTCCTTTTTGCTGCATCCCTGTCTTCCTCTGTCTTGATTTACTTGTTGTTGGTTTTCTGTTTCTTTGTTTGATTTGGTGGAAGACATAATCCCACGCTTCCTATGGAAAGGTTGTTGGGAGATTTTTAATGATTCCTCAATGTTAAAATGTCTATTTTTGTCTTGACACCCAACTAATATTTGTCTGAGCAAAACAGTCTAGATGAGAGAGAACTTCCCTGGAGGTCTGATGGCGTTTCTCCCTCGTCTTCTTA",  # noqa
        'chr2':
        "TTCAAATGAACTTCTGTAATTGAAAAATTCATTTAAGAAATTACAAAATATAGTTGAAAGCTCTAACAATAGACTAAACCAAGCAGAAGAAAGAGGTTCAGAACTTGAAGACAAGTCTCTTATGAATTAACCCAGTCAGACAAAAATAAAGAAAAAAATTTTAAAAATGAACAGAGCTTTCAAGAAGTATGAGATTATGTAAAGTAACTGAACCTATGAGTCACAGGTATTCCTGAGGAAAAAGAAAAAGTGAGAAGTTTGGAAAAACTATTTGAGGAAGTAATTGGGGAAAACCTCTTTAGTCTTGCTAGAGATTTAGACATCTAAATGAAAGAGGCTCAAAGAATGCCAGGAAGATACATTGCAAGACAGACTTCATCAAGATATGTAGTCATCAGACTATCTAAAGTCAACATGAAGGAAAAAAATTCTAAAATCAGCAAGAGAAAAGCATACAGTCATCTATAAAGGAAATCCCATCAGAATAACAATGGGCTTCTCAGCAGAAACCTTACAAGCCAGAAGAGATTGGATCTAATTTTTGGACTTCTTAAAGAAAAAAAAACCTGTCAAACACGAATGTTATGCCCTGCTAAACTAAGCATCATAAATGAAGGGGAAATAAAGTCAAGTCTTTCCTGACAAGCAAATGCTAAGATAATTCATCATCACTAAACCAGTCCTATAAGAAATGCTCAAAAGAATTGTAAAAGTCAAAATTAAAGTTCAATACTCACCATCATAAATACACACAAAAGTACAAAACTCACAGGTTTTATAAAACAATTGAGACTACAGAGCAACTAGGTAAAAAATTAACATTACAACAGGAACAAAACCTCATATATCAATATTAACTTTGAATAAAAAGGGATTAAATTCCCCCACTTAAGAGATATAGATTGGCAGAACAGATTTAAAAACATGAACTAACTATATGCTGTTTACAAGAAACTCATTAATAAAGACATGAGTTCAGGTAAAGGGGTGGAAAAAGATGTTCTACGCAAACAGAAACCAAATGAGAGAAGGAGTAGCTATACTTATATCAGATAAAGCACACTTTAAATCAACAACAGTAAAATAAAACAAAGGAGGTCATCATACAATGATAAAAAGATCAATTCAGCAAGAAGATATAACCATCCTACTAAATACATATGCACCTAACACAAGACTACCCAGATTCATAAAACAAATACTACTAGACCTAAGAGGGATGAGAAATTACCTAATTGGTACAATGTACAATATTCTGATGATGGTTACACTAAAAGCCCATACTTTACTGCTACTCAATATATCCATGTAACAAATCTGCGCTTGTACTTCTAAATCTATAAAAAAATTAAAATTTAACAAAAGTAAATAAAACACATAGCTAAAACTAAAAAAGCAAAAACAAAAACTATGCTAAGTATTGGTAAAGATGTGGGGAAAAAAGTAAACTCTCAAATATTGCTAGTGGGAGTATAAATTGTTTTCCACTTTGGAAAACAATTTGGTAATTTCGTTTTTTTTTTTTTCTTTTCTCTTTTTTTTTTTTTTTTTTTTGCATGCCAGAAAAAAATATTTACAGTAACT",  # noqa
    }

    def setup_method(self):
        self.file = pysam.FastaFile(os.path.join(BAM_DATADIR, "ex1.fa"))

    def testFetch(self):
        for id, seq in list(self.sequences.items()):
            assert self.file.fetch(id) == seq
            for x in range(0, len(seq), 10):
                assert self.file.fetch(id, x, x + 10) == seq[x:x + 10]
                # test x:end
                assert self.file.fetch(id, x) == seq[x:]
                # test 0:x
                assert self.file.fetch(id, None, x) == seq[:x]

        # unknown sequence raises IndexError
        with pytest.raises(KeyError): self.file.fetch("chr12")

    def testOutOfRangeAccess(self):
        '''test out of range access.'''
        # out of range access returns an empty string
        for contig, s in self.sequences.items():
            assert self.file.fetch(contig, len(s), len(s) + 1) == ""

    def testFetchErrors(self):
        with pytest.raises(ValueError): self.file.fetch()
        with pytest.raises(ValueError): self.file.fetch("chr1", -1, 10)
        with pytest.raises(ValueError): self.file.fetch("chr1", 20, 10)
        with pytest.raises(KeyError): self.file.fetch("chr3", 0, 100)

    def test_fetch_with_region_and_contig_raises_exception(self):
        with pytest.raises(ValueError): self.file.fetch("chr1", 10, 20, "chr1:11-20")

    def test_fetch_with_region_is_equivalent(self):
        assert self.file.fetch("chr1", 10, 20) == self.file.fetch(region="chr1:11-20")

    def testLength(self):
        assert len(self.file) == 2

    def testSequenceLengths(self):
        assert self.file.get_reference_length("chr1") == 1575
        assert self.file.get_reference_length("chr2") == 1584

    def teardown_method(self):
        self.file.close()


class TestFastaFilePathIndex:

    filename = os.path.join(BAM_DATADIR, "ex1.fa")
    data_suffix = ".fa"

    def test_raise_exception_if_index_is_missing(self):
        with pytest.raises(IOError):
            pysam.FastaFile(self.filename, filepath_index=f"garbage{self.data_suffix}.fai")

    def test_open_file_without_index_succeeds(self):
        with pysam.FastaFile(self.filename) as inf:
            assert len(inf) == 2

    def test_open_file_with_explicit_index_succeeds(self):
        with pysam.FastaFile(self.filename,
                             filepath_index=self.filename + ".fai") as inf:
            assert len(inf) == 2

    def test_open_file_with_explicit_abritrarily_named_index_succeeds(self):
        tmpfilename = get_temp_filename(self.data_suffix)
        shutil.copyfile(self.filename, tmpfilename)

        filepath_index = self.filename + ".fai"
        filepath_index_compressed = self.filename + ".gzi"
        if not os.path.exists(filepath_index_compressed):
            filepath_index_compressed = None
        with pysam.FastaFile(tmpfilename,
                             filepath_index=filepath_index,
                             filepath_index_compressed=filepath_index_compressed) as inf:
            assert len(inf) == 2

        # index should not be auto-generated
        assert not os.path.exists(tmpfilename + ".fai")
        os.unlink(tmpfilename)


class TestFastaFilePathIndexCompressed(TestFastaFilePathIndex):

    filename = os.path.join(BAM_DATADIR, "ex1.fa.gz")
    data_suffix = ".fa.gz"


class TestFastxFileFastq:

    filetype = pysam.FastxFile
    filename = "faidx_ex1.fq"
    persist = True

    def setup_method(self):
        self.file = self.filetype(os.path.join(BAM_DATADIR, self.filename),
                                  persist=self.persist)
        self.has_quality = self.filename.endswith('.fq')

    def teardown_method(self):
        self.file.close()

    def checkFirst(self, s):
        # test first entry
        assert s.sequence == "GGGAACAGGGGGGTGCACTAATGCGCTCCACGCCC"
        assert s.name == "B7_589:1:101:825:28"
        if self.has_quality:
            assert s.quality == "<<86<<;<78<<<)<;4<67<;<;<74-7;,;8,;"
            assert list(s.get_quality_array()) == [ord(x) - 33 for x in s.quality]
            assert str(s) == ("@B7_589:1:101:825:28\n"
                              "GGGAACAGGGGGGTGCACTAATGCGCTCCACGCCC\n"
                              "+\n"
                              "<<86<<;<78<<<)<;4<67<;<;<74-7;,;8,;")
        else:
            assert s.quality is None
            assert s.get_quality_array() is None
            assert str(s) == (">B7_589:1:101:825:28\n"
                              "GGGAACAGGGGGGTGCACTAATGCGCTCCACGCCC")

    def checkLast(self, s):
        assert s.sequence == "TAATTGAAAAATTCATTTAAGAAATTACAAAATAT"
        assert s.name == "EAS56_65:8:64:507:478"
        if self.has_quality:
            assert s.quality == "<<<<<;<<<<<<<<<<<<<<<;;;<<<;<<8;<;<"
            assert list(s.get_quality_array()) == [ord(x) - 33 for x in s.quality]
        else:
            assert s.quality is None
            assert s.get_quality_array() is None

    def testCounts(self):
        assert len([x for x in self.file]) == 3270

    def testMissingFile(self):
        with pytest.raises(IOError): self.filetype("nothere.fq")

    def testSequence(self):
        first = self.file.__next__()
        self.checkFirst(first)
        for last in self.file:
            pass
        self.checkLast(last)

        # test for persistence
        if self.persist:
            self.checkFirst(first)
        else:
            self.checkLast(first)

    def testManager(self):
        with self.filetype(os.path.join(BAM_DATADIR, self.filename),
                           persist=self.persist) as inf:
            first = inf.__next__()
            self.checkFirst(first)
            for last in inf:
                pass
            self.checkLast(last)

        assert inf.closed


# Test for backwards compatibility
class TestFastqFileFastq(TestFastxFileFastq):
    filetype = pysam.FastqFile


# Test for backwards compatibility
class TestFastxFileFasta(TestFastxFileFastq):
    filetype = pysam.FastqFile
    filename = "faidx_ex1.fa"


class TestFastxFileFastqStream(TestFastxFileFastq):
    persist = False


class TestFastxFileWithEmptySequence:
    """see issue 204:

    iteration over fastq file with empty sequence stops prematurely
    """

    filetype = pysam.FastxFile
    filename = "faidx_empty_seq.fq.gz"

    def testIteration(self):
        fn = os.path.join(BAM_DATADIR, self.filename)

        with gzip.open(fn) as inf:
            ref_num = len(list(inf)) / 4

        with self.filetype(fn) as f:
            l = len(list(f))
        assert ref_num == l


class TestRemoteFileFTP:
    '''test remote access.
    '''

    url = ("ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/"
           "GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa")

    def testFTPView(self):
        if not check_url(self.url):
            return

        try:
            with pysam.Fastafile(self.url) as f:
                assert len(f.fetch("chr1", 0, 1000)) == 1000
        except OSError:
            pass

    def test_sequence_lengths_are_available(self):
        if not check_url(self.url):
            return

        try:
            with pysam.Fastafile(self.url) as f:
                assert len(f.references) == 3366
                assert "chr1" in f.references
                assert f.lengths[0] == 248956422
                assert f.get_reference_length("chr1") == 248956422
        except OSError:
            pass


class TestFastqRecord:

    filetype = pysam.FastxFile
    filename = "faidx_ex1.fq"

    def setup_method(self):

        with self.filetype(os.path.join(BAM_DATADIR, self.filename), persist=True) as inf:
            self.record = next(inf)

    def test_fastx_record_sequence_can_be_modified(self):
        old_sequence = self.record.sequence
        new_record = copy.copy(self.record)
        new_sequence = "AAAC"
        new_record.set_sequence(new_sequence)
        assert str(new_record) == f">{self.record.name}\n{new_sequence}"
        assert self.record.sequence == old_sequence
        assert new_record.sequence == new_sequence

    def test_fastx_record_name_can_be_modified(self):
        old_name = self.record.name
        new_name = "new_name"
        new_record = copy.copy(self.record)
        new_record.set_name(new_name)
        assert new_record.name == new_name
        assert self.record.name == old_name

    def test_fastx_record_fail_if_name_is_None(self):
        with pytest.raises(ValueError):
            self.record.set_name(None)

    def test_fastx_record_comment_can_be_modified(self):
        old_comment = self.record.comment
        new_comment = "this is  a new comment"
        new_record = copy.copy(self.record)
        new_record.set_comment(new_comment)
        assert new_record.comment == new_comment
        assert self.record.comment == old_comment

    def test_fastx_record_comment_can_be_None(self):
        old_comment = self.record.comment
        new_comment = None
        new_record = copy.copy(self.record)
        new_record.set_comment(new_comment)
        assert new_record.comment == new_comment
        assert self.record.comment == old_comment

    def test_fastx_record_quality_can_be_modified(self):
        old_quality = self.record.quality
        new_quality = "A" * len(old_quality)
        new_record = copy.copy(self.record)
        new_record.set_sequence(self.record.sequence, new_quality)
        assert new_record.quality == new_quality
        assert self.record.quality == old_quality

    def test_fastx_record_fail_if_quality_is_wrong_length(self):
        with pytest.raises(ValueError):
            self.record.set_sequence(self.record.sequence, self.record.quality * 2)

    def test_fastx_record_can_be_created_from_scratch(self):
        fastx_record = pysam.FastxRecord()
        with pytest.raises(ValueError): str(fastx_record)
        fastx_record.set_name("name")
        with pytest.raises(ValueError): str(fastx_record)
        fastx_record.set_sequence("sequence")
        assert str(fastx_record) == ">name\nsequence"


class TestFastqProxy:

    def test_fastq_proxy_instantiation_raises_error(self):
        with pytest.raises(ValueError):
            pysam.FastqProxy()
