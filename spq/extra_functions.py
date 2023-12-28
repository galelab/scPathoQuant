__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "subprocess functions"

import subprocess
import os

def _check_subprocess_run(returncode, stderrdata, runinfo):
    if returncode == 0:
        print("STATUS: "+runinfo+" complete")
    else:
        print("WARNING: Issue with "+runinfo+" reads")
        print(stderrdata)

def _run_subprocesses(args, status, error_message, verbose=False):
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    print(status)
    stdoutdata, stderrdata = process.communicate()
    # if verbose is True:
    # print (stdoutdata)
    # print (stderrdata)
    # print (process.returncode)
    _check_subprocess_run(process.returncode, stderrdata, error_message)

