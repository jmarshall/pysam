cdef extern from "version.h":
    const char *BCFTOOLS_VERSION

def py_bcftools():
    pass

__version__ = BCFTOOLS_VERSION.decode('ascii')
