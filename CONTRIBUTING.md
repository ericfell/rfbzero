:tada: First off, thank you for taking the time to contribute! :tada:

Everyone is welcome to contribute and all types of contributions are encouraged and valued.
This includes code reviews, bug fixes, new features, examples, documentation, community participation, etc.

This document outlines the guidelines for contributing to the various aspects of the project.

# Contributing

If you find a bug in the code, a mistake in the [documentation](https://rfbzero.readthedocs.io/en/latest/index.html), or would like to see a new feature, please help us by creating [an issue in our repository](https://github.com/ericfell/rfbzero/issues), or even submit a pull request.

> And if you like the project, but just don't have time to contribute, there are other ways to support the project and show your appreciation, which we would also be very happy about:
> - Star the project
> - Tweet about it
> - Refer this project in your project's readme
> - Tell your friends/colleagues about the project



# Development Guide

## Repository Setup

1.  To work on the rfbzero.py package, you should first fork the rfbzero repository.

2.  You can then clone the fork to your computer

```bash
git clone https://github.com/<GitHubUsername>/rfbzero.py
```

3.  Make your changes and commit them to your fork.

4.  [Submit a Pull Request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests) (make sure to write an informative message so the reviewer can understand what you're adding!) via GitHub.


## Continuous Integration

`rfbzero.py` uses [GitHub Actions](https://docs.github.com/en/actions) for Continuous Integration testing. Every time you submit a pull request, a series of tests will be run to make sure the changes don’t accidentally introduce any bugs :bug:. *Your PR will not be accepted until it passes all of these tests.* You can also run tests locally to speed up the review process.

We use [Pylint](https://pylint.readthedocs.io/en/stable/) to test for [PEP 8](https://peps.python.org/pep-0008/) conformance.
To run the Pylint code analysis locally:

```
pip install pylint
pylint src/
```
:warning: if there is any output here, fix the errors and try running pylint again.


We use [mypy](https://mypy-lang.org/index.html) to type check code with type hints conforming to [PEP 484](https://peps.python.org/pep-0484/).
To run the mypy code analysis locally:

```
pip install mypy
mypy src/
```
:warning: if there is any output here, fix the errors and try running mypy again.

### Testing

`rfbzero.py` aims to have complete test coverage of our package code. If you're adding a new feature, consider writing the test first and then the code to ensure it passes. PRs which decrease code coverage will need to add tests before they can be merged.

We use [pytest](https://docs.pytest.org/en/8.2.x/) with [pytest-cov](https://github.com/pytest-dev/pytest-cov) to test code and determine code coverage.
To run tests:

```
pip install pytest pytest-cov
pytest --cov
```
:warning: if there is any output here, fix the errors and/or add tests before running pytest again.


## Acknowledgements

This CONTRIBUTING.md file was adapted from the [impedance.py GitHub repo](https://github.com/ECSHackWeek/impedance.py)

