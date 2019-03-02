#!/bin/bash

set -e;

if [ -d $TRAVIS_BUILD_DIR/build/wheelhouse ]; then
  twine upload --username pairinteraction-slave --password $PYPI_TOKEN --skip-existing --repository-url https://test.pypi.org/legacy/ $TRAVIS_BUILD_DIR/build/wheelhouse/*.whl;
fi