# SPDX-FileCopyrightText: 2024 Sebastian Weber, Henri Menke, and contributors. All rights reserved.
# SPDX-License-Identifier: GPL-3.0-or-later
from queue import Queue

# FIXME once we assume PyQt5.version >= 5.15, we can remove the except
# https://www.riverbankcomputing.com/static/Docs/PyQt5/incompatibilities.html#importing-the-sip-module
try:
    from PyQt5 import sip
except ImportError:
    import sip
from PyQt5.QtCore import QThread, pyqtSignal


class Worker(QThread):
    criticalsignal = pyqtSignal(str)

    def __init__(self, all_queues, parent=None) -> None:
        super().__init__(parent)
        self.all_queues = all_queues
        self.exiting = False

    def __del__(self) -> None:
        self.exiting = True
        if not sip.isdeleted(self):
            self.wait()

    def execute(self, stdout) -> None:
        self.stdout = stdout
        self.start()

    def run(self) -> None:
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
    def __init__(self) -> None:
        self.dataqueues = [Queue(), Queue(), Queue()]  # field1, field2, potential
        self.clear()

    def clear(self) -> None:
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

    def processOneLine(self, line) -> None:
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
            self.status_progress = (
                f"diagonalize {self.dim} x {self.dim} matrix, {self.current} of {self.total} matrices processed"
            )

        elif line[:5] == ">>OUT":
            self.current += 1
            self.status_progress = (
                f"diagonalize {self.dim} x {self.dim} matrix, {self.current} of {self.total} matrices processed"
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
