cdef extern from "version.h":
    const char *SAMTOOLS_VERSION

def py_samtools():
    pass

__version__ = SAMTOOLS_VERSION.decode('ascii')
