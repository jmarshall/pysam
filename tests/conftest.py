import http.server
import multiprocessing
import os
import socketserver
import sys
from pathlib import Path

import pytest


def _httpd(wconn, directory):
    class QuietHTTPServer(http.server.HTTPServer):
        def server_bind(self):  # Work around actions/runner-images#14568
            socketserver.TCPServer.server_bind(self)
            self.server_name, self.server_port = self.server_address[:2]

        def handle_error(self, request, client_address):
            if isinstance(sys.exc_info()[1], (BrokenPipeError, ConnectionResetError)):
                pass
            else:
                super().handle_error(request, client_address)

    class QuietRequestHandler(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=directory, **kwargs)

        def log_message(self, format, *args):
            pass

    server = QuietHTTPServer(("localhost", 0), QuietRequestHandler)
    wconn.send((server.server_name, server.server_port))
    server.serve_forever()


@pytest.fixture(scope="session")
def httpserver():
    rconn, wconn = multiprocessing.Pipe(duplex=False)
    process = multiprocessing.Process(target=_httpd, args=[wconn, Path(__file__).parent], name="httpd", daemon=True)
    process.start()

    host, port = rconn.recv()
    yield f"{host}:{port}"

    process.terminate()
    process.join()


def pytest_report_header(config):
    text = []

    if "REF_PATH" in os.environ:
        text.append("pysam: overriding REF_PATH to disable external reference lookups")
    os.environ["REF_PATH"] = ":"

    return text
