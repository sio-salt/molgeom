from collections import UserList


class _FancyIndexingList(UserList):
    def __getitem__(self, key):
        if isinstance(key, (list, tuple)):
            return _FancyIndexingList([self.data[i] for i in key])
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        if isinstance(key, (list, tuple)):
            if len(key) != len(value):
                raise ValueError("Key and value must have the same length")
            for k, v in zip(key, value):
                self.data[k] = v
        else:
            super().__setitem__(key, value)

    def __delitem__(self, key):
        if isinstance(key, (list, tuple)):
            for k in sorted(key, reverse=True):
                del self.data[k]
        else:
            super().__delitem__(key)

    def __repr__(self):
        return f"_FancyIndexingList({self.data})"
