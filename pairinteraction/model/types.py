"""Define some useful types and unions."""

from typing import Literal

from pydantic import confloat, conint, conlist, constr

PositiveFloat = confloat(ge=0)
PositiveInt = conint(ge=0)
HalfInt = confloat(multiple_of=0.5)
PositiveHalfInt = confloat(multiple_of=0.5, ge=0)

ConstituentString = Literal["atom1", "atom2", "classical_light1", "classical_light2"]
SpeciesString = constr(pattern=r"^[A-Z][a-z]?[0-9]{0,3}(_singlet|_triplet|_mqdt)?$")

Vector = conlist(float, min_length=3, max_length=3)

Symmetry = Literal[None, 1, -1]
