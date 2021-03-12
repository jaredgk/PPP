import os
import sys
import argparse
import itertools
import copy
import shutil
import logging
import subprocess


from structure import *

# Import basic structure functions
#from structure import *
from logging_module import initLogger

def run ():

    structure_err = call_structure()
    produce_structure_log(structure_err)

if __name__ == "__main__":
    initLogger()
    run()
