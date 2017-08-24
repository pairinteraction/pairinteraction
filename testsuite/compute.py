import calc.pairinteraction_complex as pi
import multiprocessing
import json
import tempfile
import os
import zmq

# Otherwise we do not find the database
os.chdir("../")

_,path_config = tempfile.mkstemp()
path_cache  = tempfile.mkdtemp()

with open(path_config, 'w') as io:
    json.dump({
        "conserveM": True,
        "dd": True,
        "deltaEPair": 1.5198298459979761e-06,
        "deltaESingle": 6.0793193839919045e-06,
        "deltaJPair": -1,
        "deltaJSingle": -1,
        "deltaLPair": -1,
        "deltaLSingle": 3,
        "deltaMPair": -1,
        "deltaMSingle": -1,
        "deltaNPair": -1,
        "deltaNSingle": 3,
        "diamagnetism": True,
        "dq": False,
        "exponent": 3,
        "invE": True,
        "invO": True,
        "j1": 0.5,
        "j2": 1.5,
        "l1": 0,
        "l2": 1,
        "m1": 0.5,
        "m2": 0.5,
        "maxBx": 0.0,
        "maxBy": 0.0,
        "maxBz": 0.0,
        "maxEx": 0.0,
        "maxEy": 0.0,
        "maxEz": 0.0,
        "maxR": 37794.52250915656,
        "minBx": 0.0,
        "minBy": 0.0,
        "minBz": 0.0,
        "minEx": 0.0,
        "minEy": 0.0,
        "minEz": 0.0,
        "minR": 377945.2250915656,
        "missingCalc": True,
        "missingWhittaker": False,
        "n1": 80,
        "n2": 79,
        "perE": True,
        "perO": True,
        "precision": 1e-12,
        "qq": False,
        "refE": False,
        "refO": False,
        "samebasis": True,
        "sametrafo": True,
        "species1": "Rb",
        "species2": "Rb",
        "steps": 10,
        "zerotheta": True
    }, io)

context = zmq.Context()
socket = context.socket(zmq.SUB)
socket.bind("tcp://*:*")
socket.setsockopt_string(zmq.SUBSCRIBE, u"")
endpoint = socket.getsockopt_string(zmq.LAST_ENDPOINT)

class Communicator:
    def __iter__(self):
        return self

    def __next__(self):
        string = socket.recv_string()

        if ">>END" in string:
            raise StopIteration
        else:
            return string

# start thread that collects the output
proc = multiprocessing.Process(
    target=pi.compute,args=(path_config, path_cache, endpoint))
proc.start()
for line in Communicator():
    print(line)

socket.close()
context.destroy()
