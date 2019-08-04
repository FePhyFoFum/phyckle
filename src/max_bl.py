#!/usr/bin/env python3

import sys
import tree_reader
import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", help="Directory containing the gene trees", required=True)
    parser.add_argument("-o","--outfile",help="Outfile for analysis",required=True)
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")
    args = parser.parse_args()

    d = args.dir
    if d[-1] != "/":
        d += "/"
    outf = open(args.outfile,"w")
    print ("gene,maxbl",file=outf)
    for i in os.listdir(d):
        #iqtree
        #if ".treefile" not in i:
        #    continue
        #raxml
        if "fastTree" not in i:
            continue
        mltree = d+i
        maxb = 0
        t = tree_reader.read_tree_file_iter(mltree).next()
        for j in t.iternodes():
            if j.length > maxb:
                maxb = j.length
        print (i.replace("RAxML_fastTreeSH_Support.","")+","+str(maxb),file=outf)
    outf.close()
