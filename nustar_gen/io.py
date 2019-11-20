import os
import warnings


def find_be_arf():
    curdir = os.path.dirname(__file__)
    tempdir = os.path.join(curdir, 'templates')
    be_arf_file = os.path.join(tempdir, 'be.arf')
    return be_arf_file
    
def find_arf_template():
    curdir = os.path.dirname(__file__)
    tempdir = os.path.join(curdir, 'templates')
    be_arf_file = os.path.join(tempdir, 'nu40101010002_srcA_sr.arf')
    return be_arf_file