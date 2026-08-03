#!/usr/bin/env python3
import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.parse_args()

params = json.loads(__import__("subprocess").check_output(["python3", "source_case.py"], cwd=__import__("os").path.dirname(__file__)))
params["m"] = 7

print(json.dumps(params))
