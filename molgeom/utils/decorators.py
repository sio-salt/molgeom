from functools import wraps
from collections.abc import Iterable


def args_to_set(func):
    """
    Decorator to normalize the arguments of a method/function to always be a set.
    """

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if (
            len(args) == 1
            and isinstance(args[0], Iterable)
            and not isinstance(args[0], (str, bytes))
        ):
            args = set(args[0])
        elif len(args) == 1 and isinstance(args[0], str):
            args = {args[0]}
        elif len(args) > 1 and any(
            isinstance(arg, (list, tuple, set, dict)) for arg in args
        ):
            raise ValueError("Multiple arguments must be passed as a single iterable")
        else:
            args = set(args)

        return func(self, args, **kwargs)

    return wrapper


def args_to_list(func):
    """
    Decorator to normalize the arguments of a method/function to always be a list.
    """

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        # # check if the first argument is an class or instance of a class
        # print(args)
        # first_arg = args[0]
        # is_class_method = inspect.isclass(first_arg)
        #
        # if is_class_method:
        #     args = args[1:]

        if (
            len(args) == 1
            and isinstance(args[0], Iterable)
            and not isinstance(args[0], (str, bytes))
        ):
            args = list(args[0])
        elif len(args) == 1 and isinstance(args[0], (str, bytes)):
            args = [args[0]]
        elif len(args) > 1 and any(
            isinstance(arg, (list, tuple, set, dict)) for arg in args
        ):
            raise ValueError("Multiple arguments must be passed as a single iterable")
        else:
            args = list(args)

        return func(self, args, **kwargs)

    return wrapper
