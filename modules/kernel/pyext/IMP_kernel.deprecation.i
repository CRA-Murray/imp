%pythoncode %{

import functools
import contextlib

def deprecated_module(version, module, help_message):
    """Mark a Python module as deprecated.
       @note The `module` argument would normally be `__name__`.
       @see [deprecation support](@ref deprecation)."""
    handle_use_deprecated(
                "Module %s is deprecated. %s" % (module, help_message))

class __deprecation_base(object):
    def __init__(self, version, help_message):
        self.version, self.help_message = version, help_message

class deprecated_object(__deprecation_base):
    """Decorator to mark a Python class as deprecated.
       @see [deprecation support](@ref deprecation)."""
    def __call__(self, obj):
        orig_init = obj.__init__
        @functools.wraps(orig_init)
        def __init__(obj, *args, **keys):
            handle_use_deprecated("Object %s is deprecated. %s"
                                  % (type(obj), self.help_message))
            orig_init(obj, *args, **keys)
        obj.__init__ = __init__
        return obj


class deprecated_method(__deprecation_base):
    """Decorator to mark a Python method as deprecated.
       @see [deprecation support](@ref deprecation)."""
    def __call__(self, obj):
        @functools.wraps(obj)
        def wrapper(cls, *args, **keys):
            handle_use_deprecated("Method %s in %s is deprecated. %s"
                               % (obj.__name__, type(cls), self.help_message))
            return obj(cls, *args, **keys)
        return wrapper


class deprecated_function(__deprecation_base):
    """Decorator to mark a Python function as deprecated.
       @see [deprecation support](@ref deprecation)."""
    def __call__(self, obj):
        @functools.wraps(obj)
        def wrapper(*args, **keys):
            handle_use_deprecated("Function %s is deprecated. %s"
                                  % (obj.__name__, self.help_message))
            return obj(*args, **keys)
        return wrapper

@contextlib.contextmanager
def allow_deprecated(allow=True):
    """Context manager to temporarily allow (or disallow) deprecated code.
       @see [deprecation support](@ref deprecation)."""
    old = get_deprecation_exceptions()
    set_deprecation_exceptions(not allow)
    yield
    set_deprecation_exceptions(old)
%}