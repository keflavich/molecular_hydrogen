"""
``--remote-data`` is provided by ``pytest-remotedata`` (an astropy plugin) when
it is installed.  Fall back to a minimal implementation so that the tests can
be run without it.
"""
import pytest

HAVE_REMOTEDATA = True
try:
    import pytest_remotedata  # noqa: F401
except ImportError:
    HAVE_REMOTEDATA = False


def pytest_addoption(parser):
    if not HAVE_REMOTEDATA:
        parser.addoption("--remote-data", action="store_true", default=False,
                         help="run tests that query remote services (HITRAN)")


def pytest_configure(config):
    config.addinivalue_line("markers",
                            "remote_data: test requires network access")


def pytest_collection_modifyitems(config, items):
    if HAVE_REMOTEDATA or config.getoption("--remote-data"):
        return
    skip = pytest.mark.skip(reason="needs --remote-data")
    for item in items:
        if "remote_data" in item.keywords:
            item.add_marker(skip)
