
import logging

class ParserVerbosity:
    SILENT = 0
    BASIC = 1
    EVERYTHING = 2

    from_string_map = {
        "silent": SILENT,
        "basic": BASIC,
        "everything": EVERYTHING
    }

def conditional_log(sender, message, level=logging.INFO, min_verbosity=ParserVerbosity.BASIC):
    """Log a message based on the verbosity level."""
    if sender.verbosity >= min_verbosity:
        sender.logger.log(level, message)