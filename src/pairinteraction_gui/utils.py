import os
from functools import wraps
from typing import Callable, Literal, Union

from PySide6.QtWidgets import QMessageBox

from pairinteraction import (
    _wrapped as pi,
    complex as pi_complex,
    real as pi_real,
)


class DatabaseMissingError(Exception):
    pass


class NoStateFoundError(Exception):
    pass


def catch_download_missing(func: Callable) -> Callable:
    @wraps(func)
    def wrapper_func(*args, **kwargs):  # noqa
        try:
            return func(*args, **kwargs)
        except RuntimeError as err:
            if "Try setting download_missing to true" in str(err):
                msg_box = QMessageBox()
                msg_box.setWindowTitle("Download missing databases?")
                msg_box.setText(str(err))
                msg_box.setInformativeText("Would you like to download the missing database?")
                msg_box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)

                if msg_box.exec() == QMessageBox.Yes:
                    set_global_database(download_missing=True)
                    return func(*args, **kwargs)
                else:
                    raise DatabaseMissingError(str(err)) from err
            raise err

    return wrapper_func


def set_global_database(
    download_missing: bool = False,
    wigner_in_memory: bool = True,
    database_dir: Union[str, "os.PathLike[str]"] = "",
) -> None:
    pi.Database._global_database = pi.Database(download_missing, wigner_in_memory, database_dir)


@catch_download_missing
def get_ket_atom(species: str, **qns: Union[float, int]) -> pi.KetAtom:
    try:
        return pi.KetAtom(species, **qns)  # type: ignore
    except ValueError as err:
        if "No state found" in str(err) or "quantum number m must be" in str(err):
            raise NoStateFoundError(str(err)) from err
        raise err


def get_basis_atom(
    ket: pi.KetAtom, *, dtype: Literal["real", "complex"] = "real", **delta_qns: Union[float, int]
) -> Union[pi_real.BasisAtom, pi_complex.BasisAtom]:
    pi = pi_real if dtype == "real" else pi_complex
    qns = {}
    for key, value in delta_qns.items():
        if value < 0:
            continue
        key = key.replace("Δ", "")
        qn = getattr(ket, key)
        qns[key] = (qn - value, qn + value)

    return pi.BasisAtom(ket.species, **qns)
