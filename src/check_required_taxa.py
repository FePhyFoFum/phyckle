#!/usr/bin/env python3

import sys
import os
import seq
import argparse


"""
this is intended to check to see if particular files have the required taxa

the file (required.txt) should look like this, each line will have a list of taxa (comma seperated). 
if one is there, it is required, if multiple, then one of the listed is required

output is
gene,bad (bad = 1 is bad)
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--required", help="File containing the required taxa.",required=True)
    parser.add_argument("-d", "--dir", help="Directory containing the gene trees", required=True)
    parser.add_argument("-o","--outfile",help="Outfile for analysis",required=True)
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")
    args = parser.parse_args()

    req = [] # list of sets
    rf = open(args.required,"r")
    for i in rf:
        spls = i.strip().split(",")
        rs = set()
        for j in spls:
            rs.add(j)
        req.append(rs)
    rf.close()

    outf = open(args.outfile,"w")

    print ("gene,required",file=outf)
    sdir = args.dir
    if sdir[-1] != "/":
        sdir+="/"
    for i in os.listdir(sdir):
        sn = i
        ss = set()
        for j in seq.read_fasta_file_iter(sdir+i):
            ss.add(j.name)
        bad = 0
        for x in req:
            if len(ss.intersection(x)) == 0:
                bad = 1
                break
        print (sn+","+str(bad),file=outf)
    outf.close()
