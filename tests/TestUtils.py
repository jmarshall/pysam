import os
import gzip
import subprocess
import time
from itertools import zip_longest

import pysam

BAM_DATADIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                           "pysam_data"))

TABIX_DATADIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                             "tabix_data"))

CBCF_DATADIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                            "cbcf_data"))

LINKDIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "linker_tests"))


def force_str(s):
    try:
        return s.decode('ascii')
    except AttributeError:
        return s


def force_bytes(s):
    try:
        return s.encode('ascii')
    except AttributeError:
        return s


def openfile(fn):
    if fn.endswith(".gz"):
        try:
            return gzip.open(fn, "rt", encoding="utf-8")
        except TypeError:
            return gzip.open(fn, "r")
    else:
        return open(fn)


def slurp_file(filename, omit_startswith=None, omit=None):
    with openfile(filename) as f:
        if omit is not None:
            return [line for line in f if not omit(line)]
        else:
            return f.readlines()


def checkBinaryEqual(filename1, filename2):
    '''return true if the two files are binary equal.
    '''
    if os.path.getsize(filename1) != os.path.getsize(filename2):
        return False

    infile1 = open(filename1, "rb")
    infile2 = open(filename2, "rb")

    def chariter(infile):
        while 1:
            c = infile.read(1)
            if c == b"":
                break
            yield c

    found = False
    for c1, c2 in zip_longest(chariter(infile1), chariter(infile2)):
        if c1 != c2:
            break
    else:
        found = True

    infile1.close()
    infile2.close()
    return found


def checkGZBinaryEqual(filename1, filename2):
    '''return true if the decompressed contents of the two files
    are binary equal.
    '''
    with gzip.open(filename1, "rb") as infile1:
        d1 = infile1.read()
        with gzip.open(filename2, "rb") as infile2:
            d2 = infile2.read()
        if d1 == d2:
            return True
    return False


def check_samtools_view_equal(
        filename1, filename2,
        without_header=False):
    '''return true if the two files are equal in their
    content through samtools view.
    '''
    # strip MD and NM tags, as not preserved in CRAM files
    args = ["-x", "MD", "-x", "NM"]
    if not without_header:
        args.append("-h")

    lines1 = pysam.samtools.view(*(args + [filename1]))
    lines2 = pysam.samtools.view(*(args + [filename2]))

    if len(lines1) != len(lines2):
        return False

    if lines1 != lines2:
        # line by line comparison
        # sort each line, as tags get rearranged between
        # BAM/CRAM
        for n, pair in enumerate(zip(lines1, lines2)):
            l1, l2 = pair
            l1 = sorted(l1[:-1].split("\t"))
            l2 = sorted(l2[:-1].split("\t"))
            if l1 != l2:
                print("mismatch in line %i" % n)
                print(l1)
                print(l2)
                return False
        else:
            return False

    return True


def dict_of_read(read, exclude=frozenset()):
    '''return a read in dictionary form, omitting excluded fields.'''
    d = {}

    # add the . for refactoring purposes.
    for x in (".query_name",
              ".query_sequence",
              ".flag",
              ".reference_id",
              ".reference_start",
              ".mapping_quality",
              ".cigartuples",
              ".next_reference_id",
              ".next_reference_start",
              ".template_length",
              ".query_length",
              ".query_qualities",
              ".bin",
              ".is_paired", ".is_proper_pair",
              ".is_unmapped", ".is_mapped",
              ".mate_is_unmapped", ".mate_is_mapped",
              ".is_reverse", ".is_forward",
              ".mate_is_reverse", ".mate_is_forward",
              ".is_read1", ".is_read2",
              ".is_secondary", ".is_qcfail",
              ".is_duplicate"):
        n = x[1:]
        if n not in exclude:
            d[n] = getattr(read, n)

    return d


def make_data_files(directory):
    if os.path.exists(os.path.join(directory, 'all.stamp')):
        return

    make = os.environ.get('MAKE', 'make')

    for attempt in range(1, 6):
        try:
            os.mkdir(os.path.join(directory, 'all.lock'), 0o700)
            break
        except FileExistsError:
            time.sleep(attempt)
            continue
    else:
        raise RuntimeError(f'Directory {directory!r} already locked: try `{make} clean` there')

    try:
        subprocess.check_output([make, '-C', directory], stderr=subprocess.STDOUT, encoding='ascii')
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f'Making test data in {directory!r} failed:\n{e.output}') from None
    finally:
        os.rmdir(os.path.join(directory, 'all.lock'))


def load_and_convert(filename, encode=True):
    '''load data from filename and convert all fields to string.

    Filename can be either plain or compressed (ending in .gz).
    '''
    data = []
    if filename.endswith(".gz"):
        with gzip.open(filename) as inf:
            for line in inf:
                line = line.decode("ascii")
                if line.startswith("#"):
                    continue
                d = line.strip().split("\t")
                data.append(d)
    else:
        with open(filename) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                d = line.strip().split("\t")
                data.append(d)

    return data


def flatten_nested_list(l):
    return [i for ll in l for i in ll]
