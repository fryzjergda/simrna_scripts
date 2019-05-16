#!/usr/bin/python

# created: 2015.08.26, Michal Boniecki
# script calculated distance restraints for SimRNA for atoms with occupancy with value different than 1.00

import sys

TOLERANCE = 0.2
WEIGHT = 1.0

if(len(sys.argv)<2):
    print >>sys.stderr, "usage: pdb_2_SimRNA_dist_restrs.py file.pdb"
    sys.exit(1)

occupancy_set = set(["1","1.0", "1.00", "1.000", "1.0000"])

input_filename = sys.argv[1]

inpfile = open(input_filename)

selected_atoms = []

curr_line = inpfile.readline()
while curr_line:

    if(curr_line[0:4]=="ATOM"):
	curr_line = curr_line.rstrip()
	curr_occupancy = curr_line[54:60].rstrip(" ").lstrip(" ")
	#print "_%s_" % curr_occupancy
	#if((curr_occupancy in occupancy_set) == True):
	#print "_%s_" % curr_occupancy
	selected_atoms.append(curr_line)
    curr_line = inpfile.readline()
inpfile.close()

for atom_line_1 in selected_atoms:
    for atom_line_2 in selected_atoms:

	if(selected_atoms.index(atom_line_1) < selected_atoms.index(atom_line_2)):

	    residue_number_1 = atom_line_1[22:26].rstrip(" ").lstrip(" ")
	    residue_number_2 = atom_line_2[22:26].rstrip(" ").lstrip(" ")

	    if(residue_number_1 != residue_number_2):

		chain_id_1 = atom_line_1[21:22]
		chain_id_2 = atom_line_2[21:22]

		atom_name_1 = atom_line_1[12:16].rstrip(" ").lstrip(" ")
		atom_name_2 = atom_line_2[12:16].rstrip(" ").lstrip(" ")

		x_1 = float(atom_line_1[30:38])
		y_1 = float(atom_line_1[38:46])
		z_1 = float(atom_line_1[46:54])

		x_2 = float(atom_line_2[30:38])
		y_2 = float(atom_line_2[38:46])
		z_2 = float(atom_line_2[46:54])

		dist = ( (x_1-x_2)**2 + (y_1-y_2)**2 + (z_1-z_2)**2 )**0.5
		
		min_dist = dist-TOLERANCE
		max_dist = dist+TOLERANCE

		print "SLOPE %s/%s/%s   %s/%s/%s  %f %f %f" % (chain_id_1, residue_number_1, atom_name_1, chain_id_2, residue_number_2, atom_name_2, min_dist, max_dist, WEIGHT)

#	    print "_%s_" % residue_number_1

