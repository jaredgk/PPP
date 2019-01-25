try:
    import ConfigParser
except ImportError:
    import configparser as ConfigParser
import argparse
import sys

is_python3 = (sys.version_info[0] == 3)



def parseOption(val):
    """Given a value read from a config file, returns appropriately formatted
    python variable.

    -If `val` == 'None', returns None variable
    -If `val` == 'True', returns True
    -If `val` == 'False', returns False
    -If `val` is a float, returns float-casted val
    -If `val` is an int and not a float, returns an int-casted val
    -If none of these criteria are met, returns the input value

    Parameters
    ----------
    val : str
        String from dict created by config parser

    Output
    ------
    val : multi
        Parsed config value for use with argparse


    """
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
    """Reads argument list to determine if there is a config file provided.

    Used to work with argparse. If `flag` is present in the input array,
    the next argument is assumed to be the filename for the config file.

    Parameters
    ----------
    arglist : list
        List of all arguments to input function, from either sys.argv or a
        custom list
    flag : str (optional) ('--conf')
        Command-line flag to look for in arglist to signify next element
        is filename for config options

    Returns
    -------
    config_filename : str (None)
        Name of config file if 'flag' is present in arglist, otherwise
        defaults to None
    is_first : bool
        Will be true if `flag` is first element of arglist. Used for method
        to allow required args to be passed from config file.

    """
    config_filename = None
    is_first = False
    if flag in arglist:
        idx = arglist.index(flag)
        config_filename = arglist[idx+1]
        if idx == 0:
            is_first = True
        return config_filename,is_first
    return None,is_first


def defaultsDictForFunction(func_name, config_name):
    """Creates a dict with default values for a given function.

    Given an input config file and a function name, this function will
    read in all values from the given section of the config file, then
    remove the variables that are used in the config file to define
    parameters for multiple functions.

    Parameters
    ----------
    func_name : str
        Name of function section of config file. (May not necessarily
        be name of function in code)
    config_name : str
        Name of config file

    Returns
    -------
    d_func : dict
        Dictionary with all attributes for the given function. Values written
        as 'None' in config file will not be added, so normal defaults
        can be used.


    """
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

def makeRequiredList(argdict, reqlist):
    """Grabs required args from argument dict created by configparser.

    Generates a list of required parameters that will be passed in lieu of
    a regular arglist or sys.argv. Must be ordered in the same way they will
    be input to argparse.

    Parameters
    ----------
    argdict : dict
        Dictionary of arguments returned from defaultsDictForFunction.
    reqlist : list
        List of required arguments, in order that they should be passed to
        argparse.

    Returns
    -------
    req_arglist : list
        List of required argument values for given function, in order for
        passing to argparse

    Exceptions
    ----------
    Argument value in argdict for required arg is None



    """
    req_arglist = []
    for arg in reqlist:
        argval = argdict[arg]
        if argval is None:
            raise Exception("Value for required argument %s cannot be None"
                            % argval)
        req_arglist.append(argdict[arg])
    return req_arglist


def checkRequired(required_args, defaults):
    for req in required_args:
        if req in defaults and req[defaults] is not None:
            raise Exception(('Required argument %s has value present '
                            ' in config file (%s) and command line' %
                            (req, defaults[req])))


def getArgsWithConfig(parser, sys_args, required_args, func_name):
    config_name, req_included = getConfigFilename(sys_args)
    if config_name is not None:
        defaults = defaultsDictForFunction(func_name, config_name)
        if not req_included:
            checkRequired(required_args, defaults)
        parser.set_defaults(**defaults)
    if not req_included:
        return parser.parse_args(sys_args)
    else:
        req_args = makeRequiredList(defaults, required_args)
        return parser.parse_args(req_args)
