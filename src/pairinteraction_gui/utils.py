# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later


class DatabaseMissingError(Exception):
    def __init__(self, err: RuntimeError) -> None:
        super().__init__(str(err))
        table = next(w for w in str(err).split(" ") if "states" in w)
        self.species = table.replace("_states", "")


class NoStateFoundError(Exception):
    def __init__(self, err: ValueError) -> None:
        super().__init__(str(err))


def get_custom_error(err: Exception) -> Exception:
    """Get a custom error message based on the type of error."""
    if isinstance(err, RuntimeError) and "Table" in str(err) and "not found" in str(err):
        return DatabaseMissingError(err)
    if isinstance(err, ValueError) and ("No state found" in str(err) or "quantum number m must be" in str(err)):
        return NoStateFoundError(err)
    return err
