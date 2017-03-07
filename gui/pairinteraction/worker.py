# Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from PyQt5 import QtCore
from queue import Queue


class Worker(QtCore.QThread):
    criticalsignal = QtCore.pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.exiting = False
        self.samebasis = False
        self.message = ""
        self.basisfile_field1 = ""
        self.basisfile_field2 = ""
        self.basisfile_potential = ""
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
        self.basisfile_field1 = ""
        self.basisfile_field2 = ""
        self.basisfile_potential = ""

        # Clear data queue
        with self.dataqueue_field1.mutex:
            self.dataqueue_field1.queue.clear()
        with self.dataqueue_field2.mutex:
            self.dataqueue_field2.queue.clear()
        with self.dataqueue_potential.mutex:
            self.dataqueue_potential.queue.clear()

        # Parse stdout
        dim = 0
        type = 0
        current = 0
        total = 0

        status_type = ""
        status_progress = ""
        status_dimension = ""

        for line in self.stdout:
            if self.exiting or not line:
                break

            elif line[:5] == u">>TYP":
                type = int(line[5:12])
                status_type = ["Field map of first atom: ",
                               "Field map of second atom: ", "Pair potential: ", "Field maps: "][type]
                status_progress = "construct matrices"

                if type == 3:
                    self.samebasis = True
                elif type == 0 or type == 1:
                    self.samebasis = False

            elif line[:5] == u">>BAS":
                basissize = int(line[5:12])
                status_progress = "construct matrices using {} basis vectors".format(
                    basissize)

            elif line[:5] == u">>STA":
                filename = line[6:-1].strip()
                if type == 0 or type == 3:
                    self.basisfile_field1 = filename
                elif type == 1:
                    self.basisfile_field2 = filename
                elif type == 2:
                    self.basisfile_potential = filename

            elif line[:5] == u">>TOT":
                total = int(line[5:12])
                current = 0

            elif line[:5] == u">>DIM":
                dim = int(line[5:12])
                status_progress = "diagonalize {} x {} matrix, {} of {} matrices processed".format(
                    dim, dim, current, total)

            elif line[:5] == u">>OUT":
                current += 1
                status_progress = "diagonalize {} x {} matrix, {} of {} matrices processed".format(
                    dim, dim, current, total)

                filenumber = int(line[5:12])
                filestep = int(line[12:19])
                blocks = int(line[19:26])
                blocknumber = int(line[26:33])
                filename = line[34:-1].strip()

                if type == 0 or type == 3:
                    self.dataqueue_field1.put(
                        [filestep, blocks, blocknumber, filename])
                elif type == 1:
                    self.dataqueue_field2.put(
                        [filestep, blocks, blocknumber, filename])
                elif type == 2:
                    self.dataqueue_potential.put(
                        [filestep, blocks, blocknumber, filename])

            elif line[:5] == u">>ERR":
                self.criticalsignal.emit(line[5:].strip())

            elif line[:5] == u">>END":
                finishedgracefully = True
                break

            else:
                print(line, end="")

            self.message = status_type + status_progress
            print("DEBUG:",self.message)

        # Clear data queue if thread has aborted
        if not finishedgracefully:
            self.clear()
