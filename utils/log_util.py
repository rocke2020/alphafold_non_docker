import logging, os


fmt = '%(asctime)s %(filename)s %(levelname)s L %(lineno)d: %(message)s'
datefmt = '%Y-%m-%d %H:%M'


def get_logger(name=None, log_file=None, log_level=logging.DEBUG):
    """ default log level DEBUG """
    logger = logging.getLogger(name)
    logging.basicConfig(format=fmt, datefmt=datefmt)
    if log_file is not None:
        log_file_folder = os.path.split(log_file)[0]
        if log_file_folder:
            os.makedirs(log_file_folder, exist_ok=True)
        fh = logging.FileHandler(log_file, 'w', encoding='utf-8')
        fh.setFormatter(logging.Formatter(fmt, datefmt))
        logger.addHandler(fh)
    logger.setLevel(log_level)
    return logger

logger = get_logger()


_nameToLevel = {
#    'CRITICAL': logging.CRITICAL,
#    'FATAL': logging.FATAL,  # FATAL = CRITICAL
#    'ERROR': logging.ERROR,
#    'WARN': logging.WARNING,
   'INFO': logging.INFO,
   'DEBUG': logging.DEBUG,
}


def get_logger_full(name=None, log_file=None, log_level=logging.DEBUG, log_level_name=''):
    """ default log level DEBUG """
    logger = logging.getLogger(name)
    logging.basicConfig(format=fmt, datefmt=datefmt)
    if log_file is not None:
        log_file_folder = os.path.split(log_file)[0]
        if log_file_folder:
            os.makedirs(log_file_folder, exist_ok=True)
        fh = logging.FileHandler(log_file, 'w', encoding='utf-8')
        fh.setFormatter(logging.Formatter(fmt, datefmt))
        logger.addHandler(fh)
    if log_level_name in _nameToLevel:
        log_level = _nameToLevel[log_level_name]
    logger.setLevel(log_level)
    return logger


def log_df_basic_info(df, comments=''):
    if comments:
        logger.info(f'comments {comments}')
    logger.info(f'df.shape {df.shape}')
    logger.info(f'df.columns {df.columns.to_list()}')
    logger.info(f'df.head()\n{df.head()}')
    logger.info(f'df.tail()\n{df.tail()}')
