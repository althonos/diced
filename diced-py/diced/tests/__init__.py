from . import test_scan


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_scan))
    return suite
