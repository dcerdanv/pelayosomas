import glob
import os
import pandas as pd
import numpy as np
from snakemake.utils import validate, min_version
#from os import listdir
#from os.path import isfile, join


#### GLOBAL scope functions ####
def get_resource(rule,resource) -> int:
	'''
	Attempt to parse config.yaml to retrieve resources available for a given
	rule. It will revert to default if a key error is found. Returns an int.
	with the allocated resources available for said rule. Ex: "threads": 1
	'''

	try:
		return config['resources'][rule][resource]
	except KeyError: # TODO: LOG THIS
		print(f'Failed to resolve resource for {rule}/{resource}: using default parameters')
		return config["resources"]['default'][resource]


#### GLOBAL PARAMETERS ####

min_version('6.2.1')

configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

OUTDIR = config['outdir']
LOGDIR = config['logdir']
DATA = config['samples']

#### Load rules ####
include: 'rules/to_long_format.smk'


def get_processed_chr():

	return [f"{OUTDIR}/split/chr21", f"{OUTDIR}/split/chr24"]



# TARGET RULES

rule all:
	input:
		get_processed_chr()
