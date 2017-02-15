#!/bin/sh

set -e


if [ "${TRAVIS_PULL_REQUEST}" = "true" ]; then
    echo "INFO: This is a PR.";
    echo "INFO: Not building docs.";
    exit 0;
fi;


if [ "${TRAVIS_BRANCH}" != "master" ]; then
    echo "INFO: We are not on the master branch.";
    echo "INFO: Not building docs.";
    #exit 0;
fi;


function gh_pages_prepare()
{
    mkdir -p "${TRAVIS_BUILD_DIR}/build/doc/doxygen/";
    cd "${TRAVIS_BUILD_DIR}/build/doc/doxygen/";

    git clone -b gh-pages https://${GH_TOKEN}@github.com/${TRAVIS_REPO_SLUG} html;

    cd "html";

    git config user.name "Travis CI";
    git config user.email "travis@travis-ci.org";

    rm -f .git/index;
    git clean -df;
}


function gh_pages_generate()
{
    cd "${TRAVIS_BUILD_DIR}/build";

    cmake -DWITH_BACKEND=No -DWITH_GUI=No -DWITH_DATABASE=No ..

    make doc
}


function gh_pages_update()
{
    cd "${TRAVIS_BUILD_DIR}/build/doc/doxygen/html";

    touch .nojekyll;

    git add --all;

    git commit -m "Doxygen build from Travis build ${TRAVIS_BUILD_NUMBER} for commit ${TRAVIS_COMMIT}"

    git push -f origin gh-pages
}


cd "${TRAVIS_BUILD_DIR}";
gh_pages_prepare;
gh_pages_generate;
gh_pages_update;
