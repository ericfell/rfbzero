#!/usr/bin/env bash
# File : run_tests.sh
# Description: Test suite driver script
set -e

# list of test case to run
tests=(
  test_redox_flow_cell.py
  test_crossover.py
  test_degradation.py
  test_experiment.py
)

# Add the module source path because we use `import rfbzero.<modulename>` in
# our test suite.  Necessary if you want to test in your local
# development environment without properly installing the package.
# this assumes we have all "test_..py" files in the 'test' folder
export PYTHONPATH="$PWD/../rfbzero" #:${PYTHONPATH}

# decide what driver to use (depending on arguments given)
if [[ $# -gt 0 && ${1} == 'coverage' ]]; then
    driver="${@} -m unittest"
elif [[ $# -gt 0 && ${1} == 'pytest' ]]; then
    driver="${@}"
elif [[ $# -gt 0 && ${1} == 'CI' ]]; then
    # Assumes the package has been installed and dependencies resolved.  This
    # would be the situation for a user.  Uses `pytest` for testing.
    shift
    unset PYTHONPATH
    driver="pytest ${@}"
else
    driver="python ${@} -m unittest"
fi

# run the tests
${driver} "${tests[@]}"