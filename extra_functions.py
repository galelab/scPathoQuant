__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "subprocess functions"

import subprocess

def _check_subprocess_run(returncode, stderrdata, runinfo):
    if returncode == 0:
        print("STATUS: "+runinfo+" complete")
    else:
        print("WARNING: Issue with "+runinfo+" reads")
        print(stderrdata)

def _run_subprocesses(args, status, error_message):
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    print("STATUS: "+status)
    stdoutdata, stderrdata = process.communicate()
    print (stdoutdata)
    _check_subprocess_run(process.returncode, stderrdata, error_message)