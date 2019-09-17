#!/bin/sh

set -e;


cat > /tmp/docker.env <<EOF
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
EOF


case "${TRAVIS_OS_NAME}" in

    # Build on GNU/Linux in Docker
    "linux")

            case "${image}" in
                "debian")
                    docker run --env-file /tmp/docker.env \
                        -v ${TRAVIS_BUILD_DIR}:/travis -w /travis \
                        --interactive --tty \
                        pairinteraction/$image \
                        /bin/bash -c "
                            set -e;
                            mkdir -p build;
                            cd build;
                            /bin/bash ../ci/prepare_documentation.sh;
                            cmake -DWITH_COVERAGE=On -DCPACK_PACKAGE_FILE_NAME=\"${package}\" ..;
                            make -k -j 2;
                            make -k -j 2 check;
                            /bin/bash ../ci/update_documentation.sh;
                        "
                    ;;

                "static-analysis")
                    docker run --env-file /tmp/docker.env \
                        -v ${TRAVIS_BUILD_DIR}:/travis -w /travis \
                        --interactive --tty \
                        pairinteraction/$image \
                        /bin/bash -c "
                            set -e;
                            /bin/bash /travis/ci/fix_style.sh;
                            mkdir -p build;
                            cd build;
                            cmake -DWITH_CLANG_TIDY=On -DWITH_GUI=Off ..;
                            make -k -j 2;
                            make -k -j 2 check;
                        "
                    ;;

                "manylinux")
                    mkdir -p build/pairinteraction_gui/pairinteraction
                    pyuic5 --output build/pairinteraction_gui/pairinteraction/plotter.py pairinteraction_gui/plotter.ui
                    docker run --env-file /tmp/docker.env \
                        -v ${TRAVIS_BUILD_DIR}:/travis -w /travis \
                        --interactive --tty \
                        pairinteraction/$image \
                        /bin/bash -c "
                            set -e;
                            mkdir -p build;
                            cd build;
                            cmake -DPYTHON_INCLUDE_DIR=\${PYTHON_INCLUDE_DIR} -DPYTHON_LIBRARY=/make/cmake/happy/ ..;
                            make -k -j 2;
                            make -k -j 2 check;
                            python setup.py bdist_wheel --python-tag py3 --plat-name manylinux2010_x86_64;
                            auditwheel repair dist/*.whl --plat manylinux2010_x86_64;
                        "
                    ;;

                *)
                    docker run --env-file /tmp/docker.env \
                        -v ${TRAVIS_BUILD_DIR}:/travis -w /travis \
                        --interactive --tty \
                        pairinteraction/$image \
                        /bin/bash -c "
                            set -e;
                            mkdir -p build;
                            cd build;
                            cmake -DCPACK_PACKAGE_FILE_NAME=\"${package}\" ..;
                            make -k -j 2;
                            make -k -j 2 check;
                            make -k -j 2 package;
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
        python setup.py bdist_wheel --python-tag py3 --plat-name macosx_10_12_x86_64;
        wheel unpack dist/*.whl -d wheelhouse/;
        chmod +x wheelhouse/*/pairinteraction/pairinteraction-* wheelhouse/*/*.data/scripts/start_pairinteraction_gui;
        cd wheelhouse/*/pairinteraction;
        python ../../../../apple/standalone.py .libs _picomplex.so _pireal.so pairinteraction-complex pairinteraction-real;
        cd ../../..;
        wheel pack wheelhouse/* -d wheelhouse/;
        ;;

esac;

if [ "${image}" = "debian" ]; then
    curl -s https://codecov.io/bash | bash -s - -X gcov
fi
