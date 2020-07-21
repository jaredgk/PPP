import os
import copy
import logging
import argparse

def local_executable (executable):

    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'bin', executable)

def confirm_executable (executable):

    '''
        Confirm if an executable exists.

        Parameters
        ----------
        executable : str
            Executable to confirm
    '''

    # Loop potental locations of executables
    for path in os.environ['PATH'].split(os.pathsep):
        
        # Join current path and executable
        executable_file = os.path.join(path, executable)
       
        # Check if the executable path exists, and if so, is an executable
        if os.path.isfile(executable_file) and os.access(executable_file, os.X_OK):

            logging.info('Calling executable: %s' % executable_file)

            # Return the path if the executable was found
            return executable_file

    # Return None if the executable was not found
    return None


def argprase_kwargs (kwarg_dict, argparse_function, bad = False):

    # Convert the kwargs
    kwargs_list = []
    for arg, value in kwarg_dict.items():
        if isinstance(value, (bool)):
            kwargs_list.append('--%s' % arg.replace('_','-'))
        elif isinstance(value, (list)) and isinstance(value[0], (list)):
            for sub_list in value: kwargs_list.extend(['--%s' % arg.replace('_','-')] + [str(sub_value) for sub_value in sub_list])
        elif isinstance(value, (list)):
            kwargs_list.extend(['--%s' % arg.replace('_','-')] + [str(sub_value) for sub_value in value])
        else:
            kwargs_list.extend(['--%s' % arg.replace('_','-'), str(value)])
    return argparse_function(kwargs_list)