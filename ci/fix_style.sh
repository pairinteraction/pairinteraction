#!/bin/bash

CLANG_FORMAT_VER=6.0
AUTOPEP8_VER=1.3.5

if ! command -v clang-format > /dev/null; then
    echo "Could not find clang-format"
    exit 1
fi

if ! command -v autopep8 > /dev/null; then
    echo "Could not find autopep8"
    exit 1
fi

if ! clang-format -version | grep -q "version ${CLANG_FORMAT_VER}"; then
    echo "Could not find clang-format ${CLANG_FORMAT_VER}!"
    echo "You have $(clang-format -version)"
    exit 1
fi

if ! autopep8 --version | grep -q "autopep8 ${AUTOPEP8_VER}"; then
    echo "Could not find autopep8 ${AUTOPEP8_VER}!"
    echo "You have $(autopep8 --version)"
    exit 1
fi

PWD=$(pwd)
SRC=$(git rev-parse --show-toplevel)

# Check if the working tree is clean
cd "${SRC}"
if ! git diff --quiet HEAD -- && [ "$1" != "-f" ]; then
    echo "Your working tree is not clean!"
    echo "Please commit your changes before continuing or use the -f option."
    exit 1
fi

# Apply clang-format
cd "${SRC}/libpairinteraction"
for file in *.h *.cpp unit_test/*.cpp; do
    echo "clang-format: $(realpath ${file})"
    clang-format -i -style=file "${file}"
done

# Apply PEP8 to Python files
cd "${SRC}/testsuite"
for file in *.py; do
    echo "PEP8: $(realpath ${file})"
    autopep8 --in-place ${file}
done
cd "${SRC}/gui/pairinteraction"
for file in *.py; do
    echo "PEP8: $(realpath ${file})"
    autopep8 --in-place ${file}
done

# Check if things have changed and print an error
cd "${SRC}"
if ! git diff --quiet HEAD -- && [ "$1" != "-f" ]; then
    echo "Formatting errors detected!"

    # Post link as a comment on GitHub
    if ! [ -z "$GH_TOKEN" ] && command -v curl > /dev/null; then
        # Upload diff
        DIFF_URL=$(git --no-pager diff | nc termbin.com 9999 | tr -d '\0')
        echo "You can download the diff from here: ${DIFF_URL}"

        # Post a comment with link to diff
        echo "Posting GitHub comment on commit ${TRAVIS_COMMIT}"
        curl --silent --show-error --output /dev/null --request POST \
             --data "{\"body\":\"Your changes have some formatting issues. You can download the diff from here: ${DIFF_URL}\n\nTo apply the diff directly use \`curl ${DIFF_URL} | git apply -\`\"}" \
             "https://api.github.com/repos/${TRAVIS_REPO_SLUG}/commits/${TRAVIS_COMMIT}/comments?access_token=${GH_TOKEN}"

        # Post failure status
        curl --silent --show-error --output /dev/null --request POST \
             --data "{\"state\": \"failed\", \"context\": \"Style check\", \"target_url\": \"${DIFF_URL}\"}" \
             "https://api.github.com/repos/${TRAVIS_REPO_SLUG}/statuses/${TRAVIS_COMMIT}?access_token=${GH_TOKEN}"
    else
        echo "Skipping GitHub comment."
    fi
else
    # Post success status
    if ! [ -z "$GH_TOKEN" ] && command -v curl > /dev/null; then
        curl --silent --show-error --output /dev/null --request POST \
             --data "{\"state\": \"success\", \"context\": \"Style check\"}" \
             "https://api.github.com/repos/${TRAVIS_REPO_SLUG}/statuses/${TRAVIS_COMMIT}?access_token=${GH_TOKEN}"
    fi
fi

cd $PWD
