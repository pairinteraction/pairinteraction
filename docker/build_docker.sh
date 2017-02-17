#!/bin/sh

set -e;


ENV_FILE="`mktemp docker.XXXXXX.env`"
SOURCE_DIR="/travis"

cat > ${ENV_FILE} <<EOF
# Travis variables
TRAVIS_COMMIT="${TRAVIS_COMMIT}"
TRAVIS_PULL_REQUEST="${TRAVIS_PULL_REQUEST}"
TRAVIS_BRANCH="${TRAVIS_BRANCH}"
TRAVIS_REPO_SLUG="${TRAVIS_REPO_SLUG}"
TRAVIS_OS_NAME="${TRAVIS_OS_NAME}"
# User variables
image="${image}"
package="${package}"
GH_TOKEN="${GH_TOKEN}"
SOURCE_DIR="${SOURCE_DIR}"
EOF


case "${TRAVIS_OS_NAME}" in

    # Build on GNU/Linux in Docker
    "linux")

            case "${image}" in
                "opensuse:13.2-w64")
                    docker run --env-file ${ENV_FILE} \
                           -v ${TRAVIS_BUILD_DIR}:${SOURCE_DIR} \
                           --interactive --tty \
                           pairinteraction/$image \
                           /bin/bash /travis/docker/build_cmake.sh win32;
                    ;;

                "debian:doxygen")
                    docker run --env-file ${ENV_FILE} \
                           -v ${TRAVIS_BUILD_DIR}:${SOURCE_DIR} \
                           --interactive --tty \
                           pairinteraction/$image \
                           /bin/bash /travis/doc/doxygen/publish_doxygen.sh
                    ;;

                *)
                    docker run --env-file ${ENV_FILE} \
                           -v ${TRAVIS_BUILD_DIR}:${SOURCE_DIR} \
                           --interactive --tty \
                           pairinteraction/$image \
                           /bin/bash /travis/docker/build_cmake.sh;
                    ;;
            esac;
            ;;

    # Build on Mac OS X directly
    "osx")

        SOURCE_DIR="."
        docker/build_cmake.sh osx;
        ;;

esac;
