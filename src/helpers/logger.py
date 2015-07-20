__author__ = 'hayssam'
import logging,os
log_dir="output/logs"
from logging.handlers import RotatingFileHandler


def setup_logger(logger, config={}):
    """
    Setup a logger from a config
    :param logger: logger to setup
    :param config: config to use
    :return:
    """
    try:
        level = config['LOGGER']['LEVEL']
    except KeyError:
        level = 'debug'
    logger.setLevel(getattr(logging, level.upper()))

    try:
        handlers = config['LOGGER']['HANDLERS']
        if isinstance(handlers, str):
            handlers = [handlers]
    except KeyError:
        handlers = ['console']

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    for h in handlers:
        hdlr = h.lower()

        if hdlr == 'console':
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(formatter)
            logger.addHandler(console_handler)

        elif hdlr == 'file':
            logfile = os.path.join(log_dir, '%s.log' % logger.name)
            try:
                max_bytes = config['LOGGER']['ROTATE_MAX_BYTES']
            except KeyError:
                max_bytes = 100000
            try:
                count = config['LOGGER']['ROTATE_COUNT']
            except KeyError:
                count = 3
            handler = RotatingFileHandler(logfile, maxBytes=max_bytes, backupCount=count)
            handler.setFormatter(formatter)
            logger.addHandler(handler)

already_set_up_logger=[]

def init_logger(name, config={}):
    global  already_set_up_logger
    logger = logging.getLogger(name='%s' % name)
    if name not in already_set_up_logger:
        setup_logger(logger, config)
        already_set_up_logger.append(name)
    return logger
