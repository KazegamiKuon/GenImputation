from contextlib import contextmanager
import os
import sys
import gzip
from itertools import groupby
import pandas as pd
import numpy as np
import subprocess
import shutil

def walk_forward(walker):
    return next(walker,False)

def get_group_file(root,ftype=None,walk_step=0,key=None):
    dir_walker = os.walk(root)
    if key is None:
        key = lambda x: x.split('_')[:-1]
    # vị trí walker đang đi tới
    walker_location = walk_forward(dir_walker)
    # số bước đã đi
    walker_step = 0
    groups = []
    while (walker_location is not False and walker_step <= walk_step):
        r, dirs, files = walker_location
        file_paths = [os.path.join(r,f) for f in files]
        file_paths = sorted(file_paths)
        groups.append(groupby(file_paths,key))
        walker_location = walk_forward(dir_walker)
        walker_step += 1
    return groups

def concat_dataframe(file_paths,header='infer',sep=','):
    '''
    input:
        file_paths: list csv file path which would be concated.
        header, sep: params as pandas.read_csv params
    output:
        return dataframe was concated
    '''
    concated_data = pd.read_csv(file_paths[0],sep=sep,header=header)
    nb_file = len(file_paths) # number csv file
    for i in range(1,nb_file):
        temp = pd.read_csv(file_paths[i],sep=sep,header=header)
        concated_data = pd.concat([concated_data,temp],ignore_index=True,join='inner')
    return concated_data

def system(command):
    subprocess.call(command, shell=True)


def system_with_stdout(command):
    proc = subprocess.Popen(
            command,
            stdin  = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
            shell  = True)
    lines = []
    for line in proc.stdout:
        lines.append(line.decode('utf-8').strip())
    return lines


def mkdir(dirname):
    if dirname.strip() != '':
        os.makedirs(dirname, exist_ok=True)


def print_error(error_message=None):
    if error_message == None:
        print(file=sys.stderr)
    else:
        print(error_message, file=sys.stderr)


@contextmanager
def reading(filename):
    root, ext = os.path.splitext(filename)
    fp = gzip.open(filename, 'rt') if ext == '.gz' else open(filename, 'rt')
    try:
        yield fp
    finally:
        fp.close()


@contextmanager
def writing(filename):
    root, ext = os.path.splitext(filename)
    fp = gzip.open(filename, 'wt') if ext == '.gz' else open(filename, 'wt')
    try:
        yield fp
    finally:
        fp.close()

def list_expression_to_dict(list_expressions):
    redict = {}
    for expr in list_expressions:
        key_values = expr.split('=')
        if len(key_values) < 2:
            redict[key_values[0]] = True
        else:
            redict[key_values[0]] = key_values[1].split(',')
    return redict

def list_to_dict(dlist: list,start_index=0):
    return dict(zip(dlist,np.arange(start_index,start_index+len(dlist))))

def line_to_items(line: str,split=' '):
    return line.rstrip().split(split)

def check_required_software(software_name):
    if os.path.exists(software_name):
        return
    if shutil.which(software_name) is not None:
        if os.path.exists(shutil.which(software_name)):
            return
    print_error()
    print_error('Error: Program "{:s}" not found in PATH'.format(
        software_name))
    print_error()
    sys.exit(0)


def check_required_file(file_path):
    if not os.path.exists(file_path):
        print_error()
        print_error('Error: File "{:s}" not found'.format(file_path))
        print_error('Please put the following file in "org_data" directory:')
        print_error()
        print_error('  ' + os.path.basename(file_path))
        print_error()
        sys.exit(0)