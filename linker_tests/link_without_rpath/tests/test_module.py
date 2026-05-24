from PysamTestModule_link_without_rpath import build_read

        
class TestModule:

    def test_pass_if_module_can_be_called(self):
        read = build_read()
        assert read.query_name == "hello"
        assert read.query_sequence == "ACGT"
