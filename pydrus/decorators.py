import os
from functools import wraps


def check_file_path(function):
    @wraps(function)
    def _check_file_path(path, *args, **kwargs):
        if not os.path.exists(path):
            raise FileNotFoundError(
                "File {} has not been found.".format(path))
        else:
            return function(path, *args, **kwargs)

    return _check_file_path
