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
                "debian")
                    docker run --env-file ${ENV_FILE} \
                        -v ${TRAVIS_BUILD_DIR}:${SOURCE_DIR} \
                        --interactive --tty \
                        pairinteraction/$image \
                        /bin/bash -c "
                            set -e;
                            cd "${SOURCE_DIR}";
                            mkdir -p build;
                            cd build;
                            /bin/bash /travis/doc/publish_documentation.sh;
                            cmake -DWITH_COVERAGE=On -DCPACK_PACKAGE_FILE_NAME="${package}" ..;
                            make -j 2;
                            make check;
                        "
                    ;;

                "ubuntu:static-analysis")
                    docker run --env-file ${ENV_FILE} \
                        -v ${TRAVIS_BUILD_DIR}:${SOURCE_DIR} \
                        --interactive --tty \
                        pairinteraction/$image \
                        /bin/bash -c "
                            set -e;
                            cd "${SOURCE_DIR}";
                            mkdir -p build;
                            cd build;
                            cmake -DWITH_CLANG_TIDY=On -DWITH_GUI=Off -DCPACK_PACKAGE_FILE_NAME="${package}" ..;
                            make -j 2;
                        "
                    ;;

                *)
                    docker run --env-file ${ENV_FILE} \
                        -v ${TRAVIS_BUILD_DIR}:${SOURCE_DIR} \
                        --interactive --tty \
                        pairinteraction/$image \
                        /bin/bash -c "
                            set -e;
                            cd "${SOURCE_DIR}";
                            mkdir -p build;
                            cd build;
                            cmake -DCPACK_PACKAGE_FILE_NAME="${package}" ..;
                            make -j 2;
                            make check;
                            make package;
                        "
                    ;;
            esac;
            ;;

    # Build on Mac OS X directly
    "osx")

        mkdir -p build;
        cd build;
        cmake -DWITH_DMG=On -DCPACK_PACKAGE_FILE_NAME="${package}" ..;
        make -j 2;
        make check;
        make package;
        make license;
        ;;

esac;

if [ "${image}" = "debian" ]; then
    cd ${TRAVIS_BUILD_DIR}
    curl -s https://codecov.io/bash | bash -
fi
