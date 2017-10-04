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

print("-- Make standalone")

executables = sys.argv[1:]
installpath = os.path.dirname(sys.argv[1])
librarypath = os.path.join(installpath, "libraries")

os.chdir(installpath)

FNULL = open(os.devnull, 'w')


def standalone(file):

    # Copy dependency
    if file not in executables:
        file_new = os.path.join(librarypath, os.path.basename(file))
        copyfile(file, file_new)
        print("Dependency {} installed.".format(os.path.basename(file)))
    else:
        file_new = file

    # Get rpaths
    out = subprocess.check_output(['otool', '-l', file]).decode('utf-8')
    p = re.compile("LC_RPATH\n.*\n.+path (.*?) \(")
    rpaths = [m.group(1) for m in p.finditer(out)]

    # Get dependencies
    o = subprocess.Popen(['otool', '-L', file], stdout=subprocess.PIPE)
    for l in o.stdout:
        l = l.decode('utf-8')

        if l[0] == '\t':

            # Get path of dependeny
            libpath = l.strip().split(' (', 1)[0]
            libpath = libpath.replace("@loader_path", os.path.dirname(file))

            if "@rpath" in libpath:
                if not rpaths:
                    print(file, rpaths, libpath)
                    raise
                for r in rpaths:
                    testpath = libpath.replace("@rpath", r)
                    testpath = testpath.replace("@loader_path", os.path.dirname(file))
                    if os.path.isfile(testpath):
                        libpath = testpath
                        break

            if not os.path.isfile(libpath):
                print(file, rpaths, libpath)
                raise

            libpath_abs = os.path.abspath(libpath)

            # If no standard library and no installed library, update path to dependency
            if not libpath_abs.startswith("/System/Library") and not libpath_abs.startswith("/usr/lib") or "libsqlite" in libpath_abs:

                # Update paths
                if os.path.basename(libpath_abs) == "Python":
                    libpath_new = "@rpath/Python"
                elif libpath_abs in executables:
                    if file not in executables:
                        libpath_new = os.path.join("@loader_path/..", os.path.basename(libpath_abs))
                    else:
                        libpath_new = os.path.join("@loader_path", os.path.basename(libpath_abs))
                else:
                    if file not in executables:
                        libpath_new = os.path.join("@loader_path", os.path.basename(libpath_abs))
                    else:
                        libpath_new = os.path.join("@loader_path", "libraries", os.path.basename(libpath_abs))

                cmd = ['install_name_tool', '-change', libpath, libpath_new, file_new]
                if subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT):
                    raise Exception("Error updating paths.")

                # Analyze the current dependency
                if os.path.basename(libpath_abs) != "Python":
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

for file in executables:
    for rpath in ["@executable_path", "@loader_path/libraries", "/usr/local/opt/python3/Frameworks/Python.framework/Versions/3.5"]: # TODO implement a function into the gui that changes the rpath of the executables according to the used python distribution and adds the pairinteraction library to the python path
        cmd = ['install_name_tool', '-add_rpath', rpath, file]
        if subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT):
            raise Exception("Error adding rpath.")
