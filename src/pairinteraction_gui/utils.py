from functools import wraps
from typing import Any, Callable, Literal, Union

from pairinteraction import (
    _wrapped as pi,
    complex as pi_complex,
    real as pi_real,
)
from pairinteraction_gui.app import Application


class DatabaseMissingError(Exception):
    def __init__(self, err: RuntimeError) -> None:
        super().__init__(str(err))
        if not self.is_database_missing_error(err):
            raise ValueError("The message must contain 'Table' and 'not found' to be a DatabaseMissingError.")
        table = next(w for w in str(err).split(" ") if "states" in w)
        self.species = table.replace("_states", "")

    @classmethod
    def is_database_missing_error(cls, err: RuntimeError) -> bool:
        return "Table" in str(err) and "not found" in str(err)


class NoStateFoundError(Exception):
    pass


def catch_download_missing(func: Callable[..., Any]) -> Callable[..., Any]:
    @wraps(func)
    def wrapper_func(*args: Any, **kwargs: Any) -> Any:
        try:
            return func(*args, **kwargs)
        except RuntimeError as err:
            if DatabaseMissingError.is_database_missing_error(err):
                Application.signals.ask_download_database.emit(DatabaseMissingError(err).species)
                return func(*args, **kwargs)
            raise err

    return wrapper_func


@catch_download_missing
def get_ket_atom(species: str, **qns: Union[float, int]) -> pi.KetAtom:
    try:
        return pi.KetAtom(species, **qns)  # type: ignore
    except ValueError as err:
        if "No state found" in str(err) or "quantum number m must be" in str(err):
            raise NoStateFoundError(str(err)) from err
        raise err


def get_basis_atom(
    ket: pi.KetAtom, *, dtype: Literal["real", "complex"] = "real", **restrict_qns: tuple[float, float]
) -> Union[pi_real.BasisAtom, pi_complex.BasisAtom]:
    if dtype == "real":
        return pi_real.BasisAtom(ket.species, **restrict_qns)  # type: ignore [arg-type]
    return pi_complex.BasisAtom(ket.species, **restrict_qns)  # type: ignore [arg-type]
