#!/usr/bin/env python
import argparse
import json
import os
import pprint
import sys

parser = argparse.ArgumentParser()
parser.add_argument ('-j','--inspect_json', action='store', required=True, dest='json', help='json output from nextflow inspect command')
parser.add_argument ('-c','--singularity_cache', action='store', required=True, dest='sinfolder', help='path where to store the singularity images')
parser.add_argument ('-u','--image_uri', action='store', default="docker://", dest='prefix', help='specify the image URI, default is <docker://>')

args = parser.parse_args()

def get_images(jsonfile, outdir):
	execDict = {}
	f = open(jsonfile, "r")
	json_content = f.read()
	# remove the footer
	firstpar = json_content.find('{')
	lastpar = json_content.rfind('}')
	json_content_clean = json_content[firstpar-1 + 1:lastpar+1]

	print(json_content_clean, file=open('output.txt', 'a'))

	# returns JSON object as a dictionary
	data = json.loads(json_content_clean)
 	# Iterating through the json list
	for process in data['processes']:
		contname = process["container"]
		newname = contname.replace('docker://', '')
		newname = newname.replace('/', '-')
		newname = outdir + "/" + newname.replace(':', '-') + ".img"
		execDict[newname] = "singularity pull --name " + newname + " " + args.prefix + contname
 
	for new_contname in execDict:
		print("Downloading to " + new_contname)
		print("Executing " + execDict[new_contname])
		os.system(execDict[new_contname])

	# Closing file
	f.close()
	return 



if __name__ == '__main__':
	jsonfile = args.json
	outdir = args.sinfolder
	if not (os.path.exists(outdir)):
		os.mkdir(outdir)
	get_images(jsonfile, outdir)
