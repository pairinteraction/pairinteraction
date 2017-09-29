import pyinteraction.picomplex as pi
import json
import tempfile
import zmq

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

pi.compute(path_config, path_cache, endpoint)

