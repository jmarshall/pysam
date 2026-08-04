cdef int samtools_invoke(int argc, char *argv[], int stdout_fd, int stderr_fd):
    samtools_set_stdout(stdout_fd)
    samtools_set_stderr(stderr_fd)
    cdef int retval = samtools_dispatch(argc, argv)
    samtools_close_stdout()
    samtools_close_stderr()
    return retval
