#!/usr/bin/env python3

import sys
import subprocess
import os
from shutil import copyfile
import re

# Literature
# https://github.com/Kitware/CMake/blob/master/Modules/BundleUtilities.cmake
# https://github.com/Stellarium/stellarium/blob/master/util/mac_app.py
# http://thecourtsofchaos.com/2013/09/16/how-to-copy-and-relink-binaries-on-osx/
# https://cmake.org/Wiki/CMake_RPATH_handling
# https://pypi.org/project/delocate/
# https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling

print("-- Make standalone")

librarypath = os.path.abspath(sys.argv[1])
executables = [os.path.abspath(executable) for executable in sys.argv[2:]]

FNULL = open(os.devnull, 'w')


def standalone(file):

    # Copy dependency
    if file not in executables:
        file_new = os.path.join(librarypath, os.path.basename(file))
        if not os.path.isfile(file_new):
            copyfile(file, file_new)
        print("Dependency {} installed.".format(os.path.basename(file)))
    else:
        file_new = file

    # Get rpaths
    out = subprocess.check_output(['otool', '-l', file]).decode('utf-8')
    p = re.compile(r"LC_RPATH\n.*\n.+path (.*?) \(")
    rpaths = [m.group(1) for m in p.finditer(out)]
    rpaths.append(os.path.dirname(file))

    # Patch rpaths
    rpaths = [rpath.replace("@loader_path", os.path.dirname(file)) for rpath in rpaths]
    rpaths = [rpath.replace("@executable_path", os.path.dirname(file)) for rpath in rpaths]

    # Get dependencies
    o = subprocess.Popen(['otool', '-L', file], stdout=subprocess.PIPE)
    for l in o.stdout:
        l = l.decode('utf-8')

        if l[0] == '\t':

            # Get path of dependeny
            libpath_original = l.strip().split(' (', 1)[0]

            # Continue if the library references itself
            if os.path.basename(libpath_original) == os.path.basename(file):
                continue

            # Patch libpath
            libpath = libpath_original.replace("@loader_path", os.path.dirname(file))
            libpath = libpath.replace("@executable_path", os.path.dirname(file))

            if "@rpath" in libpath:
                for r in rpaths:
                    testpath = libpath.replace("@rpath", r)
                    if os.path.isfile(testpath):
                        libpath = testpath
                        break
                if "@rpath" in libpath:
                    raise RuntimeError(
                        "Error in making " +
                        file +
                        " standalone. For the library " +
                        libpath_original +
                        ", there is no suitable rpath in " +
                        str(rpaths) +
                        "!")

            if not os.path.isfile(libpath):
                raise RuntimeError(
                    "Error in making " +
                    file +
                    " standalone. Library " +
                    libpath_original +
                    " not found!")

            libpath_abs = os.path.abspath(libpath)

            # If no standard library and no installed library, update path to dependency
            if not libpath_abs.startswith(
                    "/System/Library") and not libpath_abs.startswith("/usr/lib") or "libsqlite" in libpath_abs:

                # Update paths
                if libpath_abs in executables:
                    if file not in executables:
                        libpath_new = os.path.join("@loader_path/..", os.path.basename(libpath_abs))
                    else:
                        libpath_new = os.path.join("@loader_path", os.path.basename(libpath_abs))
                else:
                    if file not in executables:
                        libpath_new = os.path.join("@loader_path", os.path.basename(libpath_abs))
                    else:
                        libpath_new = os.path.join("@loader_path", "libraries", os.path.basename(libpath_abs))

                cmd = ['install_name_tool', '-change', libpath_original, libpath_new, file_new]
                if subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT):
                    raise Exception("Error in making " + file + " standalone. Error updating paths.")

                # Analyze the current dependency
                yield libpath_abs


if not os.path.exists(librarypath):
    os.makedirs(librarypath)

# http://stackoverflow.com/questions/1517614/using-otool-recursively-to-find-shared-libraries-needed-by-an-app

need = set(executables)
done = set()

while need:
    needed = set(need)
    need = set()
    for file in needed:
        need.update(standalone(file))
    done.update(needed)
    need.difference_update(done)
