import json
import os
import subprocess
import sys

case_dir = os.path.dirname(__file__)
base_case = os.path.join(case_dir, "..", "case.py")
params = json.loads(subprocess.check_output([sys.executable, base_case], cwd=case_dir))

coarse_cells = 750
params["m"] = coarse_cells - 1
params["n"] = coarse_cells - 1

print(json.dumps(params, indent=4))
