"""Defining class CustomDict.
This dictionary class keeps track of the keys that are not used.
"""
from collections.abc import MutableMapping


class CustomDict(MutableMapping):
    """Custom dictionary, that logs all the keys that are read, set or deleted and can warn about unused keys.
    We can also implement some custom behaviour.
    """

    def __init__(self, *args, unused=None, name=None, recursive=False, flattened=False, **kwargs):
        self.name = type(self).__name__ if name is None else name
        self.store = dict(*args, **kwargs)

        if flattened:

            def flatten(d):
                new = {k: v for k, v in d.items() if not isinstance(v, dict)}
                new.update(
                    {k + "." + kk: vv for k, v in d.items() if isinstance(v, dict) for kk, vv in flatten(v).items()}
                )
                return new

            self.store = flatten(self.store)

        self.unused = set(self.keys()) if unused is None else set(unused)
        if recursive:
            for key in self:
                value = self.silent_get(key)
                if isinstance(self.silent_get(key), dict):
                    self[key] = type(self)(value, name=f"{self.name}.{key}", recursive=True)

    def __getitem__(self, key):
        self.unused.discard(key)
        return self.store[key]

    def silent_get(self, key, *args):
        return self.store.get(key, *args)

    def __setitem__(self, key, value):
        self.store[key] = value
        self.unused.add(key)

    def __delitem__(self, key):
        del self.store[key]
        self.unused.discard(key)

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def __repr__(self):
        return f"{type(self).__name__}({self.store}, unused={self.unused})"

    def copy(self):
        return type(self)(self.store, unused=self.unused, name=self.name, recursive=False)

    def as_dict_recursive(self):
        return {k: v.as_dict() if isinstance(v, type(self)) else v for k, v in self.items()}

    def get_unused_recursive(self):
        unused = self.unused.copy()
        for key, v in self.store.items():
            if isinstance(v, type(self)):
                for k in v.get_unused_recursive():
                    unused.add(key + "." + k)
        return unused
