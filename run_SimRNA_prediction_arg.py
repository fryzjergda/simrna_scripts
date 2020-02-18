#!/usr/bin/python

#created, Michal Boniecki, basically for initial version of SimRNA server for RNA-Puzzles, night: 2015.08.17-18
#this script prepares and runs jobs, there will be one more for processing results
# Modified by Tomasz Wirecki 30.11.2018
# Script was adjusted to be used with SimRNA, Argparser was added instead of sys.args
# Modified by Nithin to make it compatible to SOE and SLRUM

import sys, os, shutil
import argparse

parser = argparse.ArgumentParser(prog='run_SimRNA_prediction_arg', usage='%(prog)s [options]')
parser.add_argument("-j", "--job", required=True, dest="job", type=str,
                        help="Name of jour job [required].")
parser.add_argument("-c", "--config", required=False, dest="config", default="config.dat", type=str, 
                        help="Name of configure file [default = config.dat].")
parser.add_argument("-s", "--sequence", required=False, dest="sequence", default="", type=str,
                        help="File with your sequence [optional].")
parser.add_argument("-p", "--pdb", required=False, dest="pdb", default="", type=str,
                        help="Pdb file with desired structure [optional].")
parser.add_argument("-P", "--PDB", required=False, dest="PDB", default="", type=str,          
                        help="Pdb file with desired structure - allows freezing according to occupancy [optional].")
parser.add_argument("-S", "--secstruct", required=False, dest="ss", default="", type=str,
                        help="File with secondary structure [optional].")
parser.add_argument("-E", "--replicas", required=False, dest="replicas", default=10, type=int,
                        help="Number of replicas for each run [default = 10].")                        
parser.add_argument("-n", "--runs", required=False, dest="runs", default=4, type=int,
                        help="Number of runs [default = 4].")
parser.add_argument("-r", "--restraints", required=False, dest="restraints", default="", type=str,
                        help="File with restraints [optional].")
parser.add_argument("-Q", "--queue-system",choices=["SOE", "SLRUM"],required='--grant-name' or '-G' in sys.argv, dest="queue", default="", type=str,
                        help="Type of Queue system [default = SOE]. Grant name is mandatory argument when this argument is used")
parser.add_argument("-G", "--grant-name", required='-Q' or '--queue-system' in sys.argv , dest="grant", default="", type=str,
                        help="Name of grant under which calculations are to be performed [required]. The queue system needs to be specified along with this argument")
parser.add_argument("-q", "--queue-name", required='-Q' or '--queue-system' in sys.argv , dest="partition", default="", type=str,
                        help="Name of partition under which calculations are to be performed [required]. The queue system needs to be specified along with this argument")
args = parser.parse_args()


WORKING_DIR = "WORKING_SPACE"
SIMRNA_DATA_DIR = "../../../data"
SIMRNA_BIN_PATH = "../../../bin/SimRNA "
NUMBER_OF_RUNS = args.runs

JOB_ID_NAME = args.job
SEQ_PATH = args.sequence
SS_PATH = args.ss
RE_PATH = args.restraints
E_NUM = args.replicas
CONF_PATH = args.config
PDB_PATH = args.PDB
pdb_PATH = args.pdb
QUEUE_SYS = args.queue
GRANT = args.grant
PARTITION = args.partition

# chceck input files

if pdb_PATH != "":
    if(os.path.isfile(pdb_PATH) == False): #checking pdb_PATH
        print >>sys.stderr, "pdb location/file: "+pdb_PATH+" provided by the user doesn't exist"
        sys.exit(1)
    pdb_FILENAME = pdb_PATH.split("/")[-1] #getting file name (when path is provided)

if SS_PATH != "":
    if(os.path.isfile(SS_PATH) == False): #SS_PATH
        print >>sys.stderr, "secondary structure location/file: "+SS_PATH+" provided by the user doesn't exist"
        sys.exit(2)
    SS_FILENAME = SS_PATH.split("/")[-1] #getting file name (when path is provided)
    
if PDB_PATH != "":
    if(os.path.isfile(PDB_PATH) == False): #checking pdb_PATH
        print >>sys.stderr, "PDB location/file: "+PDB_PATH+" provided by the user doesn't exist"
        sys.exit(3)
    PDB_FILENAME = PDB_PATH.split("/")[-1] #getting file name (when path is provided)
                            

if RE_PATH != "":
    if(os.path.isfile(RE_PATH) == False): #RE_PATH
        print >>sys.stderr, "restraints location/file: "+RE_PATH+" provided by the user doesn't exist"
        sys.exit(4)
    RE_FILENAME = RE_PATH.split("/")[-1] #getting file name (when path is provided)

if CONF_PATH != "config.dat":
    if(os.path.isfile(CONF_PATH) == False): #CONF_PATH
        print >>sys.stderr, "requested file: "+CONF_PATH+" is missing (in current directory)"
        sys.exit(5)
    CONF_FILENAME = CONF_PATH.split("/")[-1] #getting file name (when path is provided)
else:
    if(os.path.isfile("config.dat") == False):
        print >>sys.stderr, "requested file: config.dat is missing (in current directory)"
        print >>sys.stderr, "script termination"
        sys.exit(6)
    CONF_FILENAME = "config.dat" #getting file name (when path is provided)

if SEQ_PATH != "":           
    if(os.path.isfile(SEQ_PATH) == False): #RE_PATH
        print >>sys.stderr, "sequence location/file: "+SEQ_PATH+" provided by the user doesn't exist"
        sys.exit(4)
    SEQ_FILENAME = SEQ_PATH.split("/")[-1] #getting file name (when path is provided)
                        

if(os.path.exists("data") == True):
    os.system (" rm -rf data")
    os.symlink(SIMRNA_DATA_DIR,"data")
else:
    os.symlink(SIMRNA_DATA_DIR,"data")
'''    
    print >>sys.stderr, "requested directory (or symlink): data is missing (in current directory)"
    print >>sys.stderr, "this directory should containt energy function terms for SimRNA"
    print >>sys.stderr, "script termination"
    sys.exit(4)
'''

if(os.path.exists(WORKING_DIR) == False):
    print >>sys.stderr, "no directory "+WORKING_DIR+" which is required"
    print >>sys.stderr, "making directory: "+WORKING_DIR
    os.mkdir(WORKING_DIR)

run_dir_names = []

for curr_num in range(1,NUMBER_OF_RUNS+1):
    run_dir_names.append("run_%02d" % curr_num)

print >>sys.stderr #just empty line

os.chdir(WORKING_DIR) #entering to working_dir
print >>sys.stderr, "in "+WORKING_DIR+" making subdirectory: "+JOB_ID_NAME+"\n"
if(os.path.exists(JOB_ID_NAME) == False):
    os.mkdir(JOB_ID_NAME)
else:
    print >>sys.stderr, "requested subdirectory: "+JOB_ID_NAME+" already exists"
    print >>sys.stderr, "check it, program termination"
    sys.exit(1)

os.chdir("..")

JOB_PATH = WORKING_DIR+"/"+JOB_ID_NAME


if(SS_PATH != ""):
    print >>sys.stderr, "copying secondary structure file: "+SS_PATH+" to: "+JOB_PATH+", just to have it there (for further checking)\n"
    shutil.copy(SS_PATH, JOB_PATH)

if(RE_PATH != ""):
    print >>sys.stderr, "copying restraints file: "+RE_PATH+" to: "+JOB_PATH+", just to have it there (for further checking)\n"
    shutil.copy(RE_PATH, JOB_PATH)

if(SEQ_PATH != ""):
    print >>sys.stderr, "copying sequence file: "+SEQ_PATH+" to: "+JOB_PATH+", just to have it there (for further checking)\n"
    shutil.copy(SEQ_PATH, JOB_PATH)

if(pdb_PATH != ""):
    print >>sys.stderr, "copying pdb file: "+pdb_PATH+" to: "+JOB_PATH+", just to have it there (for further checking)\n"
    shutil.copy(pdb_PATH, JOB_PATH)

if(PDB_PATH != ""):
    print >>sys.stderr, "copying  file: "+PDB_PATH+" to: "+JOB_PATH+", just to have it there (for further checking)\n"
    shutil.copy(PDB_PATH, JOB_PATH)

if(CONF_PATH != "config.dat"):
    print >>sys.stderr, "copying configuration file: "+CONF_PATH+" to: "+JOB_PATH+", just to have it there (for further checking)\n"
    shutil.copy(CONF_PATH, JOB_PATH)
else:        
    print >>sys.stderr, "copying config.dat to "+JOB_ID_NAME+" to have it there, because parent one can be changed\n"
    shutil.copy("config.dat", JOB_PATH)

os.chdir(JOB_PATH)

print >>sys.stderr, "in "+JOB_PATH+" making directories:"
for run_dir_name in run_dir_names:

    curr_run_nr = int(run_dir_name.replace("run_","")) #geting in from current name by removing 'run_', the rest is (should be) just number of run

    print >>sys.stderr, run_dir_name

    if(os.path.exists(run_dir_name) == False):
	os.mkdir(run_dir_name)
	os.chdir(run_dir_name)
	print >>sys.stderr, "creating symlink 'data'"
	os.symlink(SIMRNA_DATA_DIR,"data")

	SimRNA_RUN_FILENAME = "r_"+JOB_ID_NAME+"_"+run_dir_name
	print >>sys.stderr, "making SimRNA run file: "+SimRNA_RUN_FILENAME
	outfile = open(SimRNA_RUN_FILENAME, "w")
	output_name = JOB_ID_NAME+"_"+run_dir_name
	
	command = SIMRNA_BIN_PATH
	
	if CONF_PATH != "config.dat":
	    command += "-c ../" +CONF_FILENAME + " "
        else:
            command += "-c ../config.dat "
            
        if pdb_PATH != "" and PDB_PATH == "" and SEQ_PATH == "":
            command += "-p ../" +pdb_FILENAME+ " "
        elif pdb_PATH == "" and PDB_PATH != "" and SEQ_PATH == "":
            command += "-P ../" +PDB_FILENAME+ " "
        elif PDB_PATH == "" and pdb_PATH == "" and SEQ_PATH != "":
            command += "-s ../" +SEQ_FILENAME+ " "
        else:
            print >>sys.stderr, "You must provide EXACTLY one of the following files: -s, -p, -P"
            sys.exit(1)
        
        if SS_PATH != "":
            command += "-S ../" +SS_FILENAME+ " "
        
        if RE_PATH != "":
            command += "-r ../" +RE_FILENAME+ " "
        
        if E_NUM != "10":
            command += "-E " +str(E_NUM)+ " "
        else:
            command += "-E 10 "
        
        command += "-R "+str(curr_run_nr)+" -o "+output_name+" >& "+SimRNA_RUN_FILENAME+".out"            
            
        #print command
	
	if QUEUE_SYS == "SLRUM":
            	print >>outfile,"#!/bin/bash -l"
            	print >>outfile,"#SBATCH -J "+output_name
            	print >>outfile,"#SBATCH --output \""+SimRNA_RUN_FILENAME+".stdout\""
            	print >>outfile,"#SBATCH --error \""+SimRNA_RUN_FILENAME+".stderr\""
            	print >>outfile,"#SBATCH -A \""+GRANT+"\""
            	print >>outfile,"#SBATCH -p "+PARTITION
            	print >>outfile,"#SBATCH --mem 1600"
            	print >>outfile,"#SBATCH -N 1"
            	print >>outfile,"#SBATCH --tasks-per-node 1"
            	print >>outfile,"#SBATCH --cpus-per-task 10"
            	print >>outfile,"#SBATCH --time 72:00:00"

	print >> outfile, command
	
	outfile.close()
	os.system("chmod +x "+SimRNA_RUN_FILENAME)
	os.chdir("..")
    else:
	print >>sys.stderr, "directory: "+run_dir_name+" already exists, script termination"
	sys.exit(5)

print >>sys.stderr #just empty line

print >>sys.stderr, "submitting jobs:"
for run_dir_name in run_dir_names:
    if(os.path.exists(run_dir_name) == True):
    
	os.chdir(run_dir_name)

	SimRNA_RUN_FILENAME = "r_"+JOB_ID_NAME+"_"+run_dir_name #should be same as before, in previous loop
	#exit(1)
	if QUEUE_SYS == "SOE":
            	os.system("qsub -cwd -q all.q -l h_vmem=300M -l mem_free=500M -pe mpi "+str(E_NUM)+" "+SimRNA_RUN_FILENAME)  
	if QUEUE_SYS == "SLRUM":
            	os.system("sbatch -n1 --account="+GRANT+" "+SimRNA_RUN_FILENAME)
	os.chdir("..")
    else:
	print >>sys.stderr, "directory: "+run_dir_name+" doesn't exist, script termination"
	sys.exit(6)
