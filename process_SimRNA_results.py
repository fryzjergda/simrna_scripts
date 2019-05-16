#!/usr/bin/python

#created, Michal Boniecki, for automatic processing results from SimRNA runs, 2015.12.14
#this script prepares processes output data from SimRNA runs

import sys, os, shutil
from glob import glob

if(len(sys.argv) < 2):
    print >>sys.stderr, "usage: process_SimRNA_results.py job_id_name"
    print >>sys.stderr, "note: job_id_name should be legal string, to be used to set up working directory for this job"
    print >>sys.stderr, "note: job_id_name should be the same, as used in the script that launches simulations"
    sys.exit(1)

WORKING_DIR = "WORKING_SPACE"
RUN_MASK = "run_??"
PROCESSING_DIR = "processing_results"

JOB_ID_NAME = sys.argv[1]

JOB_PATH = WORKING_DIR+"/"+JOB_ID_NAME
ALL_TRAFL_filename = JOB_ID_NAME+"_ALL.trafl"

PDB_REFERENCE_FILE_FOR_TRAFL_CONVERSION = JOB_ID_NAME+"_run_01_01-000001.pdb"
PDB_REFERENCE_PATH_FOR_TRAFL_CONVERSION = "../run_01/"+PDB_REFERENCE_FILE_FOR_TRAFL_CONVERSION

ALL_TRAFL_low_size_thrs = 1024 # 1kb for the time being

FRACTION_LOWEST_ENERGY_FRAMES_TO_CLUSTER = 0.01

OUTPUT_PDBS_DIR = "output_PDBS"


if(os.path.exists(JOB_PATH) == False):
    print >>sys.stderr, "expected path: "+JOB_PATH+" doesn't exist"
    print >>sys.stderr, "this path should exist for given job_id_name: "+JOB_ID_NAME
    sys.exit(2)

dir_list = glob(JOB_PATH+"/"+RUN_MASK)
dir_list.sort()
#print dir_list

n_runs = len(dir_list)
if(n_runs == 0):
    print >>sys.stderr, "inside directory: "+JOB_PATH+" there are no expected directories named: "+RUN_MASK+" where ?? is numbering field: 01, 02, 03 ..."
    sys.exit(3)

print >>sys.stderr, "number of run directories detected: "+str(n_runs)

os.chdir(JOB_PATH)
print >>sys.stderr, "making directory: "+JOB_PATH+"/"+PROCESSING_DIR
if(os.path.exists(PROCESSING_DIR) == False):
    os.mkdir(PROCESSING_DIR)
else:
    print >>sys.stderr, "requested subdirectory: "+PROCESSING_DIR+" already exists"
    print >>sys.stderr, "check it, maybe delete it, ... program termination"
    sys.exit(4)

os.chdir(PROCESSING_DIR)
command = "cat ../"+RUN_MASK+"/*.trafl > "+ALL_TRAFL_filename
print >>sys.stderr, "being in "+JOB_PATH+"/"+PROCESSING_DIR+"  running command:"
print >>sys.stderr, command
os.system(command)

if(os.path.isfile(ALL_TRAFL_filename) == False):
    print >>sys.stderr, "expected (from previous step) file: "+ALL_TRAFL_filename+" doesn't exist"
    sys.exit(5)

file_size = os.path.getsize(ALL_TRAFL_filename)
if(file_size < ALL_TRAFL_low_size_thrs):
    print >>sys.stderr, "file: "+ALL_TRAFL_filename+" is too small: "+str(file_size)+" bytes"
    print >>sys.stderr, "it seems the file contains no data, something went wrong before ..."
    sys.exit(6)

# some tests if ALL_TRAFL_filename is correct (sometimes there are problems, when there is no disk space during SimRNA running)
# if ALL_TRAFL_filename is not correct, it should be repaired here

inpfile = open(ALL_TRAFL_filename)
first_line = inpfile.readline()
second_line = inpfile.readline().rstrip()
inpfile.close()

#assuming that second line in file ALL_TRAFL_filename is first line containing coordinated, thus is possible to calculate the size of system (seq length)
#by dividing of number of items by 15 (3 coordinates x,y,z and 5 atoms per nucleotide)

coords_list = second_line.split()
seq_length = len(coords_list) / 15

if(seq_length < 4):
    print >>sys.stderr, "it seems that seq_lenght detected from file: "+ALL_TRAFL_filename+" is too low"
    print >>sys.stderr, "something went wrong"
    sys.exit(7)

print >>sys.stderr, "clustering ... assuming:"
print >>sys.stderr, "--- fraction of lowest energy frames to clustering: "+str(FRACTION_LOWEST_ENERGY_FRAMES_TO_CLUSTER)
rmsd_thrs = 0.1*float(seq_length)
rmsd_thrs_str = "%.1f" % rmsd_thrs
print >>sys.stderr, "--- rmsd thrs for clustering 0.1*seq_lenght which is: "+rmsd_thrs_str

command = "../../../bin/clustering "+ALL_TRAFL_filename+" "+str(FRACTION_LOWEST_ENERGY_FRAMES_TO_CLUSTER)+" "+rmsd_thrs_str+" > clustering.log 2>&1"
print >>sys.stderr, command

os.system(command)

clust_1_2_3_names = glob("*clust0[1-3].trafl")
clust_1_2_3_names.sort()

n_clusts = len(clust_1_2_3_names)
print >>sys.stderr, "number of clusters to process: "+str(n_clusts)

print >>sys.stderr, "extracting pdbs, reconstructing all atom representation"
print >>sys.stderr, "creating symlink 'data'"
os.symlink("../../../data","data")

if(os.path.exists(PDB_REFERENCE_PATH_FOR_TRAFL_CONVERSION) == False):
    print >>sys.stderr, "expected file at location: "+PDB_REFERENCE_PATH_FOR_TRAFL_CONVERSION+" doesn't exist"
    sys.exit(8)

for curr_clust_name in clust_1_2_3_names:

    clust_reconstr_log_name = curr_clust_name.replace(".trafl",".log")
    command = "../../../bin/SimRNA_trafl2pdbs "+PDB_REFERENCE_PATH_FOR_TRAFL_CONVERSION+" "+curr_clust_name+" 1 AA > "+clust_reconstr_log_name+" 2>&1"
    print >>sys.stderr, command

    os.system(command)

os.chdir("..")

print >>sys.stderr, "making directory: "+OUTPUT_PDBS_DIR
if(os.path.exists(OUTPUT_PDBS_DIR) == False):
    os.mkdir(OUTPUT_PDBS_DIR)
else:
    print >>sys.stderr, "directory already exists"

os.chdir(PROCESSING_DIR)

pdbs_list = glob("*.pdb")
#pdbs_AA_list = glob("*_AA.pdb")
ss_detected_list = glob("*.ss_detected")
if(len(pdbs_list) < 0):
    print >>sys.stderr, "inside directory: "+PROCESSING_DIR+" there is no pdb files, something when wrong in previous step"
    sys.exit(9)
else:
    pdbs_list.sort()
    print >>sys.stderr, "detected pdb files in: "+PROCESSING_DIR+":"
    for curr_pdb_name in pdbs_list:
	print >>sys.stderr, curr_pdb_name

print >>sys.stderr, "copying pdb and ss_detected files to: "+OUTPUT_PDBS_DIR+" just to store them there"

for curr_pdb_name in pdbs_list:
    shutil.copy(curr_pdb_name, "../"+OUTPUT_PDBS_DIR)

for curr_ss_detected_name in ss_detected_list:
    shutil.copy(curr_ss_detected_name, "../"+OUTPUT_PDBS_DIR)

os.chdir("..")

print >>sys.stderr, "DONE :-)"
