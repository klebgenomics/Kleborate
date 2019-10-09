# Kleborate tests

Kleborate comes with a few automated tests to help with development and spotting bugs. Most users don't need to worry about these, but you're welcome to run the tests if you want!

To run the tests, first navigate to Kleborate's root directory (i.e the directory you made when cloning it from GitHub which contains `kleborate-runner.py`).

If you have [pytest](https://docs.pytest.org/en/latest/) installed, you can use it to run the tests:
```
cd Kleborate
python3 -m pytest
```

Otherwise, the tests can be run using Python's [built-in unit testing framework](https://docs.python.org/3/library/unittest.html):
```
cd Kleborate
python3 -m unittest
```
