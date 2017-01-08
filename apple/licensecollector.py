#!/usr/bin/env python3

import glob
import sys
import os

path_in = sys.argv[1]
path_out = sys.argv[2]

with open(path_out, "w") as fw:
    with open(os.path.join(path_in, "LICENSE"), "r") as fr:
        fw.write(fr.read())

    notes = """-----------------------------------------------------------------------------

    3RD PARTY LICENSES

    Note that not all files in this software package belong to the pairinteraction
    project. For 3rd party files, the individual licenses apply. The licenses
    are contained in pairinteraction.app/Contents/Resources/licenses.

    Boost Software License - Version 1.0 - 17 August 2003
    https://www.boost.org/LICENSE_1_0.txt

    GNU GENERAL PUBLIC LICENSE - Version 3 - 29 June 2007
    https://www.gnu.org/licenses/gpl-3.0.html

    GNU LESSER GENERAL PUBLIC LICENSE - Version 3 - 29 June 2007
    https://www.gnu.org/licenses/lgpl-3.0.html

    The MIT License
    https://opensource.org/licenses/MIT

    Mozilla Public License - Version 2.0
    https://www.mozilla.org/en-US/MPL/2.0/
"""

    fw.write(notes)
