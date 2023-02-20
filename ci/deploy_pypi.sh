#!/bin/bash

set -e;

if [ -d $TRAVIS_BUILD_DIR/build/wheelhouse ]; then
  case "$1" in
    "real")
      twine upload --username pairinteraction-slave --password $PYPI_TOKEN --skip-existing $TRAVIS_BUILD_DIR/build/wheelhouse/*.whl;
      ;;
    "test")
      twine upload --username pairinteraction-slave --password $PYPI_TOKEN --skip-existing --repository-url https://test.pypi.org/legacy/ $TRAVIS_BUILD_DIR/build/wheelhouse/*.whl;
      ;;
  esac;
fi
