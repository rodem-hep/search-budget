def pytest_configure(config):
    config.addinivalue_line("markers", "slow: rebuilds outputs and compares them byte for byte")
