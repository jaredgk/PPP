import sys
import logging



def formatLogger():
    fmt_def = "%(asctime)s - %(name)s - %(levelname)s: %(message)s"
    fmtr = logging.Formatter(fmt=fmt_def)
    s_handler = logging.StreamHandler()
    s_handler.setFormatter(fmtr)
    rootlogger = logging.getLogger()
    rootlogger.addHandler(s_handler)

def switchLogger(name):
    filename = name+'.log'
    fmt_def = "%(asctime)s - %(name)s - %(levelname)s: %(message)s"
    fmtr = logging.Formatter(fmt=fmt_def)
    f_handler = logging.FileHandler(filename)
    l = logging.getLogger()
    f_handler.setFormatter(fmtr)
    l.addHandler(f_handler)
