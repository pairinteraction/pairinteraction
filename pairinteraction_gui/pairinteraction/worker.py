# Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
#
# This file is part of the pairinteraction GUI.
#
# The pairinteraction GUI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The pairinteraction GUI is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the pairinteraction GUI. If not, see <http://www.gnu.org/licenses/>.
from queue import Queue

from PyQt5 import QtCore


class Worker(QtCore.QThread):
    criticalsignal = QtCore.pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.exiting = False
        self.samebasis = False
        self.message = ""
        self.basisfiles_field1 = []
        self.basisfiles_field2 = []
        self.basisfiles_potential = []
        self.dataqueue_field1 = Queue()
        self.dataqueue_field2 = Queue()
        self.dataqueue_potential = Queue()

    def __del__(self):
        self.exiting = True
        self.wait()

    def execute(self, stdout):
        self.stdout = stdout
        self.start()

    def clear(self):
        with self.dataqueue_field1.mutex:
            self.dataqueue_field1.queue.clear()
        with self.dataqueue_field2.mutex:
            self.dataqueue_field2.queue.clear()
        with self.dataqueue_potential.mutex:
            self.dataqueue_potential.queue.clear()

    def run(self):
        finishedgracefully = False

        self.message = ""

        # Clear filenames
        self.basisfiles_field1 = []
        self.basisfiles_field2 = []
        self.basisfiles_potential = []

        # Clear data queue
        with self.dataqueue_field1.mutex:
            self.dataqueue_field1.queue.clear()
        with self.dataqueue_field2.mutex:
            self.dataqueue_field2.queue.clear()
        with self.dataqueue_potential.mutex:
            self.dataqueue_potential.queue.clear()

        # Parse stdout
        dim = 0
        _type = 0
        current = 0
        total = 0

        status_type = ""
        status_progress = ""

        for line in iter(self.stdout.readline, b""):
            if isinstance(line, str):
                line = line.encode("utf-8")

            if self.exiting or not line:
                break

            elif line[:5] == b">>TYP":
                _type = int(line[5:12].decode("utf-8"))
                status_type = [
                    "Field map of first atom: ",
                    "Field map of second atom: ",
                    "Pair potential: ",
                    "Field maps: ",
                ][_type]
                status_progress = "construct matrices"

                if _type == 3:
                    self.samebasis = True
                elif _type == 0 or _type == 1:
                    self.samebasis = False

            elif line[:5] == b">>BAS":
                basissize = int(line[5:12].decode("utf-8"))
                status_progress = f"construct matrices using {basissize} basis vectors"

            elif line[:5] == b">>STA":
                filename = line[6:-1].decode("utf-8").strip()
                if _type == 0 or _type == 3:
                    self.basisfiles_field1.append(filename)
                elif _type == 1:
                    self.basisfiles_field2.append(filename)
                elif _type == 2:
                    self.basisfiles_potential.append(filename)

            elif line[:5] == b">>TOT":
                total = int(line[5:12].decode("utf-8"))
                current = 0

            elif line[:5] == b">>DIM":
                dim = int(line[5:12])
                status_progress = "diagonalize {} x {} matrix, {} of {} matrices processed".format(
                    dim, dim, current, total
                )

            elif line[:5] == b">>OUT":
                current += 1
                status_progress = "diagonalize {} x {} matrix, {} of {} matrices processed".format(
                    dim, dim, current, total
                )

                filestep = int(line[12:19].decode("utf-8"))
                blocks = int(line[19:26].decode("utf-8"))
                blocknumber = int(line[26:33].decode("utf-8"))
                filename = line[34:-1].decode("utf-8").strip()

                if _type == 0 or _type == 3:
                    self.dataqueue_field1.put([filestep, blocks, blocknumber, filename])
                elif _type == 1:
                    self.dataqueue_field2.put([filestep, blocks, blocknumber, filename])
                elif _type == 2:
                    self.dataqueue_potential.put([filestep, blocks, blocknumber, filename])

            elif line[:5] == b">>ERR":
                self.criticalsignal.emit(line[5:].decode("utf-8").strip())

            elif line[:5] == b">>END":
                finishedgracefully = True
                break

            else:
                print(line.decode("utf-8"), end="")

            self.message = status_type + status_progress

        # Clear data queue if thread has aborted
        if not finishedgracefully:
            self.clear()
        self.stdout.close()
