cdef int bcftools_invoke(int argc, char *argv[], int stdout_fd, int stderr_fd):
    bcftools_set_stdout(stdout_fd)
    bcftools_set_stderr(stderr_fd)
    cdef int retval = bcftools_dispatch(argc, argv)
    bcftools_close_stdout()
    bcftools_close_stderr()
    return retval
