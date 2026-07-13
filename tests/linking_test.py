"""test linking against pysam.
"""

import os
import pytest
import subprocess
import pysam

from TestUtils import LINKDIR


if not os.environ.get("PYSAM_LINKING_TESTS"):
    pytest.skip("enable linking tests by setting PYSAM_LINKING_TESTS environment variable", allow_module_level=True)


def check_import(statement):
    try:
        output = subprocess.check_output(
            statement, stderr=subprocess.STDOUT, shell=True)
    except subprocess.CalledProcessError as exc:
        if b"ImportError" in exc.output:
            raise ImportError(
                "module could not be imported: {}".format(str(exc.output)))
        else:
            raise


def check_pass(statement):
    try:
        output = subprocess.check_output(
            statement, stderr=subprocess.STDOUT, shell=True)
    except subprocess.CalledProcessError as exc:
        raise ValueError("{}: {}".format(exc, exc.output))
    if b"FAILED" in output:
        raise ValueError("module tests failed")
    return True


class TestLinking:
    package_name = "link_with_rpath"

    def setup_method(self):
        self.workdir = os.path.join(LINKDIR, self.package_name)
        self.testdir = os.path.join(self.workdir, "tests")

    def test_package_can_be_installed(self):
        subprocess.check_output(f"cd {self.workdir} && rm -rf build && python setup.py install", shell=True)


class TestLinkWithRpath(TestLinking):
    package_name = "link_with_rpath"

    def test_package_tests_pass(self):
        assert check_pass(f"cd {self.testdir} && python test_module.py")


class TestLinkWithoutRpath(TestLinking):
    package_name = "link_without_rpath"

    def test_package_tests_fail_on_import(self):
        with pytest.raises(ImportError):
            check_import(f"cd {self.testdir} && python test_module.py")

    def test_package_tests_pass_if_ld_library_path_set(self):
        pysam_libraries = pysam.get_libraries()
        pysam_libdirs, pysam_libs = zip(
            *[os.path.split(x) for x in pysam_libraries])
        pysam_libdir = pysam_libdirs[0]

        assert check_pass(f"export LD_LIBRARY_PATH={pysam_libdir}:$PATH && cd {self.testdir} && python test_module.py")
