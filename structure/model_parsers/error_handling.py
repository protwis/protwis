import traceback

def log_or_raise(logger, message, exception_type=Exception, action="log", parent_exception=None):
    if action in ["log", "log_then_raise"]:
        logger.error(message)

    if action in ["log_with_trace", "log_with_trace_then_raise"]:
        logger.error(message, exc_info=True)

    if action == "raise" or action == "log_then_raise":
        if parent_exception is not None:
            raise exception_type(message) from parent_exception
        else:
            raise exception_type(message)
        traceback.print_exc() #Prints a full stack trace preserving the parent exception location

    if action not in ["log", "raise", "log_then_raise", "log_with_trace", "log_with_trace_then_raise"]:
        raise ValueError(f"Invalid action '{action}' specified. Must be 'log', 'raise', or 'log_then_raise'.")