from pairinteraction import backend, gui, model, preprocessing, simulation
from pairinteraction.backend._backend import VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH
from pairinteraction.units import ureg

__all__ = [
    "gui",
    "backend",
    "model",
    "preprocessing",
    "simulation",
    "ureg",
]

__version__ = f"{VERSION_MAJOR}.{VERSION_MINOR}.{VERSION_PATCH}"
