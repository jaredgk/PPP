import sys
import logging



def individualFunctionLogger():
    fmt_def = "%(asctime)s - %(name)s - %(levelname)s: %(message)s"
    fmtr = logging.Formatter(fmt=fmt_def)
    s_handler = logging.StreamHandler()
    s_handler.setFormatter(fmtr)
    s_handler.setLevel(logging.WARNING)
    f_handler = logging.FileHandler('root.log',mode="w")
    f_handler.setFormatter(fmtr)
    rootlogger = logging.getLogger()
    rootlogger.setLevel('INFO')
    rootlogger.addHandler(s_handler)
    rootlogger.addHandler(f_handler)
    def exp_handler(etype,val,tb):
        logging.error("%s" % (val), exc_info=(etype,val,tb))
    #exp_handler.append(sys.__excepthook__)
    sys.excepthook = exp_handler

def pipeSwitchLogger(name):
    filename = name+'.log'
    fmt_def = "%(asctime)s - %(name)s - %(levelname)s: %(message)s"
    fmtr = logging.Formatter(fmt=fmt_def)
    f_handler = logging.FileHandler(filename)
    l = logging.getLogger()
    f_handler.setFormatter(fmtr)
    l.addHandler(f_handler)
