#!/usr/bin/env bash
# File       : check_coverage.sh
# Description: Coverage wrapper around test suite driver script
set -e

tool='coverage'
if [[ $# -gt 0 ]]; then
    # optional argument to use different tool to check coverage
    tool="${1}"; shift
fi

if [[ ${tool} == 'coverage' ]]; then
    # run the tests (generates coverage data to build report)
    ./run_tests.sh coverage run --source=. "${@}"

    # build the coverage report on stdout
    coverage report -m
elif [[ ${tool} == 'pytest' ]]; then
    # generate coverage reports with pytest in one go
    # ./run_tests.sh pytest --cov=../rfbzero --cov-report=term-missing "${@}"
    # when not including playground .py file
    ./run_tests.sh pytest --cov-config=../.coveragerc --cov=../rfbzero --cov-report=term-missing "${@}"

else
    # error: write to stderr
    >&2 echo "Error: unknown tool '${tool}'"
    exit 1
fi