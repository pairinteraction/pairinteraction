#!/bin/sh

set -e;


export ENV_FILE="`mktemp docker.XXXXXX.env`";
export SOURCE_DIR="/travis";


cat > ${ENV_FILE} <<EOF
# Travis variables
TRAVIS_COMMIT=${TRAVIS_COMMIT}
TRAVIS_PULL_REQUEST=${TRAVIS_PULL_REQUEST}
TRAVIS_BRANCH=${TRAVIS_BRANCH}
TRAVIS_REPO_SLUG=${TRAVIS_REPO_SLUG}
TRAVIS_OS_NAME=${TRAVIS_OS_NAME}
# User variables
image=${image}
package=${package}
GH_TOKEN=${GH_TOKEN}
SOURCE_DIR=${SOURCE_DIR}
EOF


case "${TRAVIS_OS_NAME}" in

    # Build on GNU/Linux in Docker
    "linux")

            case "${image}" in
                "debian:doxygen")
                    docker run --env-file ${ENV_FILE} \
                           -v ${TRAVIS_BUILD_DIR}:${SOURCE_DIR} \
                           --interactive --tty \
                           hmenke/$image \
                           /bin/bash /travis/doc/publish_documentation.sh
                    ;;

                "debian")
                    docker run --env-file ${ENV_FILE} \
                           -v ${TRAVIS_BUILD_DIR}:${SOURCE_DIR} \
                           --interactive --tty \
                           hmenke/$image \
                           /bin/bash /travis/docker/build_cmake.sh gcov;
                    ;;

                *)
                    docker run --env-file ${ENV_FILE} \
                           -v ${TRAVIS_BUILD_DIR}:${SOURCE_DIR} \
                           --interactive --tty \
                           hmenke/$image \
                           /bin/bash /travis/docker/build_cmake.sh;
                    ;;
            esac;
            ;;

    # Build on Mac OS X directly
    "osx")

        export SOURCE_DIR=".";
        docker/build_cmake.sh osx;
        ;;

esac;

if [ "${image}" = "debian" ]; then
    cd ${TRAVIS_BUILD_DIR}
    curl -s https://codecov.io/bash | bash -
fi
