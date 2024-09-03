# Kleborate tests

Kleborate comes with a few automated tests to help with development and spotting bugs. To run them, you'll need the `pytest` and `pytest-mock` packages installed.

To run the tests, first navigate to Kleborate's root directory (i.e the directory you made when cloning it from GitHub which contains `kleborate-runner.py`, execute this command:
```
python3 -m pytest
```

Or if you have [Coverage.py](https://coverage.readthedocs.io) installed, you can run the tests through it to get some code coverage stats:
```
coverage run --source . -m pytest && coverage report -m
```
