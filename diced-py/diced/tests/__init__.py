from . import test_scan, test_doctest


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_scan))
    suite.addTests(loader.loadTestsFromModule(test_doctest))
    return suite
