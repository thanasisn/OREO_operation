# -*- coding: utf-8 -*-
"""
General functions for controlling and logging

@author: thanasisn
"""

import os
import sys
from datetime import datetime
from dotmap   import DotMap
import yaml

def goodbye(logfile, tic, scriptname, quiet=False):
    """
    Log script execution to a central file
    """
    out  = datetime.now().strftime("%F %T")    + " "
    out += os.getlogin() + "@" + os.uname()[1] + " "
    out += os.path.normpath(scriptname)        + " "
    out += str(round((datetime.now() - tic).total_seconds() / 60.0, 2)) + " mins"
    if not quiet:
        print('\n' + out + '\n')
    with open(logfile, "a", encoding="utf-8") as runlog:
        runlog.write(out + '\n')


def size_to_human(num, suffix="B"):
    """
    Convert numbers to human readable format
    """
    for unit in ("", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"):
        if abs(num) < 1024.0:
            return f"{num:3.1f} {unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f} Yi{suffix}"


def true_if_file_exist(file, quiet=False):
    """
    Skip current iteration if file already exist
    """
    if os.path.isfile(file):
        if not quiet:
            print("Skip existing file:",
                  os.path.basename(file) + ",",
                  datetime.fromtimestamp((os.path.getmtime(file))).strftime("%F %T") + ",",
                  size_to_human(os.path.getsize(file)))
        return True
    return False


def output_needs_update(filein, fileout, quiet=False, minsize=5):
    """
    Check if we have to produce a new file if the source
    file is older or the new file is empty.

    Parameters
    ----------
    filein : string
        Source data file.

    fileout : string
        Target file.

    Returns
    -------
    Boolean
    """
    #  can not resolve missing input file
    if not os.path.isfile(filein):
        sys.exit(f"Input file {filein} does not exist!")
    #  output file is almost empty
    if os.path.getsize(fileout) < minsize:
        if not quiet:
            print(f"Update: {fileout} is only {os.path.getsize(fileout)} bytes")
        return True
    #  input is newer
    if os.path.getmtime(filein) > os.path.getmtime(filein):
        if not quiet:
            print(f"Update: {fileout} is older than {filein}")
        return True
    #  need to update outfile
    print(f"Skip: {fileout}")
    return False


def get_configs(file):
    """
    Read configuration profile file.
    This should provide all option to run the project.

    Parameters
    ----------
    file : string
        Path to a yaml file.

    Returns
    -------
    A DotMap object.
    """
    if not os.path.isfile(file):
        sys.exit(f"Missing config file: {file}")

    print("\nOpening config file:", file, "\n")
    with open(file, 'r', encoding="utf-8") as file:
        configs = yaml.safe_load(file)
        # Convert dictionary to use dot notation
        return DotMap(configs)
