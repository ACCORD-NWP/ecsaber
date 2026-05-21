#!/usr/bin/env python3
# This script creates equivalent yaml files for a use with OOPS-ECMWF

import argparse
import datetime
import os
import yaml
import copy
import sys
import collections.abc

# Correct date formatting
def correct_date(d):
    if isinstance(d, dict):
        # Loop over dictionary items
        for k,v in d.items():
            d[k] = correct_date(v)
    elif isinstance(d, list):
        # Loop over list items
        i = 0
        for v in d:
            d[i] = correct_date(v)
            i += 1

    if isinstance(d, datetime.datetime):
        # Replace with string
        d = d.strftime("%Y-%m-%dT%H:%M:%SZ")
    return d

# Replace key name recursively
def rename_key(d, key, repl):
    if isinstance(d, dict):
        # Loop over dictionary items
        newd = {}
        for k,v in d.items():
            if k == key:
                newd[repl] = v
            else:
                newd[k] = rename_key(v, key, repl)
    elif isinstance(d, list):
        # Loop over list items
        newd = []
        for v in d:
            newd.append(rename_key(v, key, repl))
    else:
        newd = d
    return newd

# Find pattern in value and replace it
def find_and_replace(d, value, repl):
    if isinstance(d, dict):
        # Loop over dictionary items
        newd = {}
        for k,v in d.items():
             newd[k] = find_and_replace(v, value, repl)
    elif isinstance(d, list):
        # Loop over list items
        newd = []
        for v in d:
            newd.append(find_and_replace(v, value, repl))
    elif isinstance(d, str):
        newd = d.replace(value, repl)
    else:
        newd = d
    return newd

# Switch hybrid and ensemble covariance models to SABER covariance model
def to_saber(config, ensembleToSaber):
    # Key: "covariance model" => "covariance"
    if "covariance model" in config: 
        config["covariance"] = config["covariance model"]
        config.pop("covariance model")

    if "covariance" in config:
        # hybrid => SABER hybrid
        if config["covariance"] == "hybrid":
              config["covariance"] = "SABER"
              config["covariance type"] = "hybrid"

        # ensemble => SABER ensemble
        if config["covariance"] == "ensemble" and ensembleToSaber:
            config["covariance"] = "SABER"
            config["covariance type"] = "ensemble"
            if "members" in config:
                config["ensemble"] = {}
                config["ensemble"]["members"] =  config["members"]
                config.pop("members")
            elif "members from template" in config:
                config["ensemble"] = {}
                config["ensemble"]["members from template"] =  config["members from template"]
                config.pop("members from template")

    if "covariance type" in config:
        # Recursive call for SABER hybrid
        if config["covariance type"] == "hybrid":
            components = []
            for cmp in config["components"]:
                newCmp = {}
                newCmp["covariance"] = to_saber(cmp["covariance"], True)
                newCmp["weight"] = cmp["weight"]
                components.append(newCmp)
            config["components"] = components

    return config

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("inputYaml", help="Yaml file name")
parser.add_argument("outputYaml", help="Yaml file name")
args = parser.parse_args()

# Read yaml file
print("--  - Updating yaml: " + args.inputYaml)
with open(args.inputYaml, "r") as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

# Correct date formatting recursively
correct_date(config)

# Rename "state variables" into "variables"
config = rename_key(config, "state variables", "variables")

# Background
if "background" in config:
    background = {}
    if "states" in config["background"]:
      background["state"] = config["background"]["states"]
    else:
      background["state"] = [config["background"]]
    config["background"] = background
    date = config["background"]["state"][0]["date"]
    variables = config["background"]["state"][0]["variables"]

# Background error
if "background error" in config:
    ensembleToSaber = ("members from template" in config["background error"])
    config["background error"] = to_saber(config["background error"], ensembleToSaber)
    if config["background error"]["covariance"] == "hybrid":
       config["background error"]["covariance"] = "SABER"
       config["background error"]["covariance type"] = "hybrid"

# Add model
config["model"] = {}
config["model"]["tstep"] = "PT6H"

# Output as Increment4D
for key in ["output dirac", "output perturbations", "output variance"]:
  if key in config:
    if "states" in config[key]:
      config[key]["increment"] = config[key]["states"]
      config[key].pop("states")

# Output as State4D
for key in ["output states"]:
  if key in config:
    if "states" in config[key]:
      config[key]["state"] = config[key]["state"]
      config[key].pop("states")

# Write yaml file
with open(args.outputYaml, "w") as file:
    output = yaml.dump(config, file, sort_keys=False)
