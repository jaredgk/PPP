try:
    import ConfigParser
except:
    import configparser as ConfigParser
import argparse
import sys

is_python3 = (sys.version_info[0] == 3)



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
        idx = arglist.index(flag)
        return arglist[idx+1]
    return None


def defaultsDictForFunction(func_name, config_name):
    if is_python3:
        config = ConfigParser.ConfigParser()
    else:
        config = ConfigParser.SafeConfigParser()
    config.read(config_name)
    d_vars = dict(config.items("DEFAULT"))
    d_hold = dict(config.items(func_name))
    d_func = {}
    for k in d_hold.keys():
        if k not in d_vars:
            val = parseOption(d_hold[k])
            if val is not None:
                d_func[k] = val
    return d_func
