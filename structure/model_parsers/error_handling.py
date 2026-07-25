def log_or_raise(logger, message, exception_type=Exception, action="log", parent_exception=None):
    if action == "log" or action == "log_then_raise":
        logger.error(message)
    if action == "raise" or action == "log_then_raise":
        raise exception_type(message) from parent_exception
    if action not in ["log", "raise", "log_then_raise"]:
        raise ValueError(f"Invalid action '{action}' specified. Must be 'log', 'raise', or 'log_then_raise'.")
