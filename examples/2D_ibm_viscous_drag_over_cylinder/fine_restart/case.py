import json
import os
import subprocess
import sys

case_dir = os.path.dirname(__file__)
base_case = os.path.join(case_dir, "..", "case.py")
params = json.loads(subprocess.check_output([sys.executable, base_case], cwd=case_dir))

params["old_grid"] = "T"
params["old_ic"] = "T"
params["t_step_old"] = 0
params["n_start"] = 0
params["num_patches"] = 0

for key in list(params):
    if key.startswith("patch_icpp("):
        del params[key]

print(json.dumps(params, indent=4))
