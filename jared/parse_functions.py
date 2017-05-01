import ConfigParser
import argparse


def getArgsFromParser(parser):
    for k in vars(parser)
    

def parseOption(val):
    if val == 'None':
        return None
    if val == 'True':
        return True
    if val == 'False':
        return False
    try:
        a = int(val)
        b = float(val)
        if a == b:
            return int(val)
        return float(val)
    except:
        pass
    return val

def getConfigFilename(arglist, flag='--conf'):
    if flag in arglist:
        idx = argslist.index(flag)
        return argslist[idx+1]
    return None
        
    
def defaultsDictForFunction(func_name, config_name):
    config = ConfigParser.SafeConfigParser()
    config.read(config_name)
    d_vars = dict(config.items("DEFAULT"))
    d_hold = dict(config.items(func_name))
    d_func = {}
    for k in d_vars.keys():
        if k not in d_func:
            val = parseOption(d_hold[k])
            if val is not None:
                d_func[k] = val
    return d_func
            