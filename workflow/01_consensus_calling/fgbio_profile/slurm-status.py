#!/usr/bin/env python3
import subprocess
import sys

jobid = sys.argv[1]

try:
    out = subprocess.check_output(
        ["squeue", "-h", "-j", jobid, "-o", "%T"]
    ).decode().strip()

    if out == "":
        print("success")
    elif out in ("RUNNING", "PENDING", "CONFIGURING", "COMPLETING"):
        print("running")
    else:
        print("failed")

except Exception:
    print("failed")
