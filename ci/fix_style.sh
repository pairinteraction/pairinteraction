#!/usr/bin/env bash

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

if ! command -v curl > /dev/null; then
    echo "Could not find curl"
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

PREVDIR=$(pwd)
SRC=$(git rev-parse --show-toplevel)

# Check if the working tree is clean
cd "${SRC}"
if ! git diff --quiet HEAD -- && [ "$1" != "-f" ]; then
    echo "Your working tree is not clean!"
    echo "Please commit your changes before continuing or use the -f option."
    exit 1
fi

# Apply clang-format
find "${SRC}/pairinteraction" -name "*.h" -or -name "*.cpp" -print0 | xargs -0 -P 0 -I '{}' -t clang-format -i -style=file "{}"

# Apply PEP8 to Python files
# FIXME: Currently ignoring E225 and E226 to prevent changing @LIBNAME@ to @ LIBNAME @
find "${SRC}/testsuite" -name "*.py" -print0 | xargs -0 -P 0 -I '{}' -t autopep8 --max-line-length 120 --aggressive --ignore E225,E226 --in-place "{}"
find "${SRC}/pairinteraction_gui/pairinteraction" -name "*.py" -print0 | xargs -0 -P 0 -I '{}' -t autopep8 --max-line-length 120 --aggressive --in-place "{}"

# Check if things have changed and print an error
cd "${SRC}"
if ! git diff --quiet HEAD -- && [ "$1" != "-f" ]; then
    echo "Formatting errors detected!"

    # Post link as a comment on GitHub
    if [ -n "$GH_TOKEN" ]; then
        # Upload diff
        DIFF_URL=$(git --no-pager diff | nc termbin.com 9999 | tr -d '\0')
        echo "You can download the diff from here: ${DIFF_URL}"

        # Post a comment with link to diff
        echo "Posting GitHub comment on commit ${TRAVIS_COMMIT}"
        curl --silent --show-error --output /dev/null --request POST \
             --header "Authorization: token ${GH_TOKEN}" \
             --data "{\"body\":\"Your changes have some formatting issues. You can download the diff from here: ${DIFF_URL}\n\nTo apply the diff directly use \`curl ${DIFF_URL} | git apply -\`\"}" \
             "https://api.github.com/repos/${TRAVIS_REPO_SLUG}/commits/${TRAVIS_COMMIT}/comments"

        # Post failure status
        curl --silent --show-error --output /dev/null --request POST \
             --header "Authorization: token ${GH_TOKEN}" \
             --data "{\"state\": \"failure\", \"context\": \"Style check\", \"target_url\": \"${DIFF_URL}\"}" \
             "https://api.github.com/repos/${TRAVIS_REPO_SLUG}/statuses/${TRAVIS_COMMIT}"
    else
        echo "I can't comment on GitHub because GH_TOKEN is missing so I'm crashing the build."
        echo "Here is what's wrong:"
        git --no-pager diff
        exit 1
    fi
else
    # Post success status
    if [ -n "$GH_TOKEN" ]; then
        curl --silent --show-error --output /dev/null --request POST \
             --header "Authorization: token ${GH_TOKEN}" \
             --data "{\"state\": \"success\", \"context\": \"Style check\"}" \
             "https://api.github.com/repos/${TRAVIS_REPO_SLUG}/statuses/${TRAVIS_COMMIT}"
    fi
fi

cd "$PREVDIR"
