"""Defining Custom exceptions for this package."""


class QnNotFoundError(Exception):
    """Exception raised when a qn is not found in a list."""

    def __init__(self, qn, useAll=False, msg=None):
        self.qn = qn
        self.useAll = useAll
        self.msg = msg

    def __str__(self):
        if self.msg:
            return self.msg
        listName = "basisQnumbers" if not self.useAll else "allQnumbers"
        return f"Quantum number {self.qn} not found in list {listName}."


class _CppDeleted:
    """Simple class to distinguish between None and deleted objects (see deleteCppObjects)."""

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(cls, cls).__new__(cls)
        return cls._instance

    def __repr__(self):
        return "This Cpp object was deleted, you can not use this object anymore."


CppDeleted = _CppDeleted()  # singleton but with a nice repr


class CppObjectAlreadyDeleted(Exception):
    """Exception raised when a object, that was already deleted is called."""

    def __init__(self, objectname, msg=None):
        self.objectname = objectname
        self.msg = msg

    def __str__(self):
        if self.msg:
            return self.msg
        return f"The cpp object {self.objectname} was already deleted and cannot be called anymore."
