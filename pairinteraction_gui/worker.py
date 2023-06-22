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

from PyQt5.QtCore import QThread, pyqtSignal


class Worker(QThread):
    criticalsignal = pyqtSignal(str)

    def __init__(self, all_queues, parent=None):
        super().__init__(parent)
        self.all_queues = all_queues
        self.exiting = False

    def __del__(self):
        self.exiting = True
        self.wait()

    def execute(self, stdout):
        self.stdout = stdout
        self.start()

    def run(self):
        self.exiting = False

        for line in iter(self.stdout.readline, b""):
            if self.exiting or not line or self.all_queues.finishedgracefully:
                break
            self.all_queues.processOneLine(line)

        # Clear data queue if thread has aborted
        if not self.all_queues.finishedgracefully:
            self.all_queues.clear()
        if self.stdout is not None:
            self.stdout.close()


class AllQueues:
    def __init__(self):
        self.dataqueues = [Queue(), Queue(), Queue()]  # field1, field2, potential
        self.clear()

    def clear(self):
        self.basisfiles = [[], [], []]  # field1, field2, potential
        self.dataqueues = [Queue(), Queue(), Queue()]  # field1, field2, potential

        # Parse stdout
        self.dim = 0
        self._type = 0
        self.current = 0
        self.total = 0

        self.message = ""
        self.status_type = ""
        self.status_progress = ""

        self.finishedgracefully = False

    def processOneLine(self, line):
        if isinstance(line, bytes):
            line = line.decode("utf-8")

        if line[:5] == ">>TYP":
            self._type = int(line[5:12])
            self.status_type = [
                "Field map of first atom: ",
                "Field map of second atom: ",
                "Pair potential: ",
                "Field maps: ",
            ][self._type]
            self.status_progress = "construct matrices"

        elif line[:5] == ">>BAS":
            basissize = int(line[5:12])
            self.status_progress = f"construct matrices using {basissize} basis vectors"

        elif line[:5] == ">>STA":
            filename = line[6:-1].strip()
            self.basisfiles[self._type % 3].append(filename)

        elif line[:5] == ">>TOT":
            self.total = int(line[5:12])
            self.current = 0

        elif line[:5] == ">>DIM":
            self.dim = int(line[5:12])
            self.status_progress = "diagonalize {} x {} matrix, {} of {} matrices processed".format(
                self.dim, self.dim, self.current, self.total
            )

        elif line[:5] == ">>OUT":
            self.current += 1
            self.status_progress = "diagonalize {} x {} matrix, {} of {} matrices processed".format(
                self.dim, self.dim, self.current, self.total
            )

            filestep = int(line[12:19])
            blocks = int(line[19:26])
            blocknumber = int(line[26:33])
            filename = line[34:-1].strip()
            self.dataqueues[self._type % 3].put([filestep, blocks, blocknumber, filename])

        elif line[:5] == ">>ERR":
            self.criticalsignal.emit(line[5:].strip())

        elif line[:5] == ">>END":
            self.finishedgracefully = True

        else:
            print(line, end="")

        self.message = self.status_type + self.status_progress
