#!/usr/bin/env python3

import sys
import os
import subprocess
import argparse

"""
this will use bp to check for whether the outgroup is monophyletic

it makes some assumptions about files. the ml files are just X.fasta.treefile
"""


# SH, 
BPCMD = "bp -t MLTREEFILE -c OUTGROUP -scut 80"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-g", "--outgroupf", help="File containing the tree with outgroup designated.",required=True)
    parser.add_argument("-d", "--dir", help="Directory containing the gene trees", required=True)
    parser.add_argument("-s", "--filenamepattern",help="is there a particular common part of the name?",
                        default="",required =False )
    parser.add_argument("-e","--exclude",help="exclude files with this string",default="",required=False)
    parser.add_argument("-o","--outfile",help="Outfile for analysis",required=True)

    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    args = parser.parse_args()
    
    outgroupf = args.outgroupf
    print ("outgroup file: "+outgroupf,file=sys.stderr)

    gened = args.dir
    if gened[-1] != "/":
        gened += "/"
    
    if args.filenamepattern != "":
        print("checking all files in '"+gened+"' with '"+args.filenamepattern+"' in the filename",file=sys.stderr)
    else:
        print("checking all files in '"+gened+"'",file=sys.stderr)

    outf = open(args.outfile,"w")

    for i in os.listdir(gened):
        if args.filenamepattern in i:
            if len(args.exclude) > 0:
                if args.exclude in i:
                    continue
            gf = i
            cmd = BPCMD.replace("MLTREEFILE",gened+i).replace("OUTGROUP",outgroupf)
            process = subprocess.Popen(cmd.split(" "),stdout = subprocess.PIPE,stderr=subprocess.PIPE)
            out,err = process.communicate()
            out = str(out)
            count = 0
            start = False
            for j in out.split("\n"):
                if len(j.strip()) > 10:
                    count += 1
            print(gf+","+str(count),file=outf)
    outf.close()