import traceback

handling_modes = ["log", "raise", "log_then_raise", "log_with_trace", "log_with_trace_then_raise"]

def log_or_raise(logger, message, exception_type=Exception, action="log", parent_exception=None):
    """
    Log a message and/or raise an exception based on the specified action.

    Parameters
    ----------
    logger : logging.Logger
        The logger instance to use for logging messages.
    message : str
        The message to log or include in the raised exception.
    exception_type : Exception, optional
        The type of exception to raise. Default is Exception.
    action : str, optional
        The action to take. Must be one of "log", "raise", "log_then_raise", "log_with_trace", or "log_with_trace_then_raise". Default is "log".
    parent_exception : Exception, optional
        The original exception to include in the raised exception. Default is None.
    """
    if action in ["log", "log_then_raise"]:
        logger.error(message)

    if action in ["log_with_trace", "log_with_trace_then_raise"]:
        logger.error(message, exc_info=True)

    if action in ["raise", "log_then_raise", "log_with_trace_then_raise"]:
        if parent_exception is not None:
            traceback.print_exception(type(parent_exception), parent_exception, parent_exception.__traceback__)
            raise exception_type(message) from parent_exception
        else:
            traceback.print_exc()
            raise exception_type(message)

    if action not in handling_modes:
        raise ValueError(f"Invalid action '{action}' specified. Must be one of { ', '.join(handling_modes)}.")