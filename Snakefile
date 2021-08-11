import glob
import os
import pandas as pd
import numpy as np
from snakemake.utils import validate, min_version
from os import listdir
from os.path import isfile, join
from pathlib import Path   


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

def get_params(rule,param) -> int:
	'''
	Attempt to parse config.yaml to retrieve parameters available for a given
	rule. It will crash otherwise.
	''' 
	try:
		return config['parameters'][rule][param]
	except KeyError: # TODO: LOG THIS
		print(f'Failed to resolve parameter for {rule}/{param}: Exiting...')
		sys.exit(1)


#### GLOBAL PARAMETERS ####

min_version('6.2.1')

configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

OUTDIR = config['outdir']
LOGDIR = config['logdir']
DATA = config['samples']
CHR = config['chr']

#### Load rules ####
include: 'rules/to_long_format.smk'


#def get_processed_chr():

	# For rule remove_header
	#return [f"{OUTDIR}/chr21.ready_noheader.vcf", f"{OUTDIR}/chr24.ready_noheader.vcf"]

	# For rule split_chr
	#return [f"{OUTDIR}/split/chr21", f"{OUTDIR}/split/chr24"]
	#return f"{OUTDIR}/split/chr24"
	
	# For rule paste_header
	#all_input = expand(f"{OUTDIR}/headed/chr21/{{part}}", part = [f for f in listdir("results/split/chr21") \
	#			if isfile(join("results/split/chr21", f)) and not f.startswith('.')])
	#all_input = expand(f"{OUTDIR}/headed/chr24/{{pedazo}}", pedazo = [f for f in listdir("results/split_aux/chr24") \
	#			if isfile(join("results/split_aux/chr24", f)) and not f.startswith('.')])
	#print(all_input)
	#return all_input

	# For rule make_long
	#all_input = expand(f"{OUTDIR}/long_part/chr24/{{part}}", part = [f for f in listdir("results/headed/chr24") \
	#				if isfile(join("results/headed/chr24", f)) and not f.startswith('.')])

	#parts = glob_wildcards(f"{OUTDIR}/header/chr24/chr24.part-{{part}}").part
	#all_input = expand(f"{OUTDIR}/long_part/chr24/chr24.part-{{partes}}", partes = parts)
	#print(all_input)
	#return all_input




########################################


#def pair_name_to_infiles():
#	# get all *.clean.vcf files recursively under DIR_A
#	vcf_path = Path('results/split/chr24').glob('**/chr24.part-*')
#
#	inputt = []
#
#	for f in vcf_path:
#		inputt.append(str(f))
#	print(inputt)
#
#	return inputt


# using function written in python code, map vcf name to their infile path
# = pair_name_to_infiles()

#rule all:
#	input:
#		pair_name_to_infiles()
        #expand('DIR_A/dir3/{vcf_name}.output.vcf', vcf_name=vcf_infiles_dict.keys())



###################################


# TARGET RULES

#rule all:
#	input:
#		get_processed_chr()

rule all:
	input:
		f"{OUTDIR}/final/{CHR}.tsv"