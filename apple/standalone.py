#!/usr/bin/env python3

import sys
import subprocess
import os
from shutil import copyfile

# Literature
# https://github.com/Kitware/CMake/blob/master/Modules/BundleUtilities.cmake
# https://github.com/Stellarium/stellarium/blob/master/util/mac_app.py
# http://thecourtsofchaos.com/2013/09/16/how-to-copy-and-relink-binaries-on-osx/
# https://cmake.org/Wiki/CMake_RPATH_handling

executables = sys.argv[1:]
installpath = os.path.dirname(sys.argv[1])
librarypath = os.path.join(installpath, "libraries")

FNULL = open(os.devnull, 'w')


def standalone(file):
    # Copy library
    if file not in executables:
        file_new = os.path.join(librarypath, os.path.basename(file))
        copyfile(file, file_new)
    else:
        file_new = file

    # Get dependencies
    o = subprocess.Popen(['otool', '-L', file], stdout=subprocess.PIPE)
    for l in o.stdout:
        l = l.decode('utf-8')
        if l[0] == '\t':
            libpath = l.strip().split(' (', 1)[0]
            libpath = libpath.replace("@loader_path", os.path.dirname(file))
            libpath = libpath.replace("@rpath", "/usr/local/opt/llvm/lib")
            if not libpath.startswith("/usr/lib") or "libsqlite" in libpath:

                # Update paths
                if file not in executables:
                    libpath_new = os.path.join("@loader_path", os.path.basename(libpath))
                else:
                    libpath_new = os.path.join("@executable_path", "libraries", os.path.basename(libpath))

                cmd = ['install_name_tool', '-change', libpath, libpath_new, file_new]
                if subprocess.call(cmd, stdout=FNULL, stderr=subprocess.STDOUT):
                    raise Exception("Error updating paths.")

                yield libpath


print("-- Make standalone")

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
