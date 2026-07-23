def log_or_raise(logger, message, exception_type=Exception, action="log", parent_exception=None):
    if action == "log":
        logger.error(message)
    else:
        raise exception_type(message) from parent_exception