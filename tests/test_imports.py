# validate that the package works and has been imported correctly


def test_imports():
    import snakealtpromoter
    import snakealtpromoter.workflows

    assert snakealtpromoter is not None
    assert snakealtpromoter.workflows is not None
