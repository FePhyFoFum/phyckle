#!/usr/bin/env python3

import sys
import os
import argparse

"""
this is meant to take the constraint file
tree
-tree
-tree

and output the csv files with whether the mltrees match the constraints
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--constraintsf", help="File containing the constraints designated.",required=True)
    parser.add_argument("-s", "--check_outgroup_mono_script_location", help="Where is the check_outgroup_mono.py script?",default="check_outgroup_mono.py",required=False)
    parser.add_argument("-d", "--dir", help="Directory containing the gene trees", required=True)

    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    args = parser.parse_args()

    cosl = args.check_outgroup_mono_script_location

    mldir = args.dir
    tf = open(args.constraintsf,"r")
    count = 0
    for i in tf:
        if i[0] == "-":
            i = i[1:].strip()
        else:
            i = i.strip()
        of = open("constraints"+str(count),"w")
        of.write(i+';')
        of.close()
        cmd = cosl+" -g constraints"+str(count)+" -d "+mldir+" -o constraints"+str(count)+".csv"
        print (cmd,file=sys.stderr)
        os.system(cmd)
        count += 1
    tf.close()
