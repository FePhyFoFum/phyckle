#!/usr/bin/env python3

import argparse
import sys
import os
import os.path

cmd = "iqtree -nt 3 -s SEQ -alrt 1000 -sp SEQ.parts -pre NAME -m GTR+G -redo"
cmd = "IQTREE -nt THREADS -s SEQ -alrt 1000 -pre NAME -m GTR+G -redo"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--iqtree", help="Where is the iqtree bin?",default="iqtree",required=False)
    parser.add_argument("-t", "--threads",help="How many threads for raxml?",default="4",required=False)
    parser.add_argument("-d", "--dir",help="Which directory?",required=True)
    parser.add_argument("-s", "--filenamepattern",help="is there a particular common part of the name?",
                        default="",required =False )

    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    args = parser.parse_args()

    std = args.dir
    isn = args.filenamepattern

    if std[-1] != "/":
        std += "/"

    cmd = cmd.replace("IQTREE",args.iqtree).replace("THREADS",args.threads)

    for i in os.listdir(std):
        if len(isn) > 0:
            if isn not in i:
                continue
        x = cmd.replace("SEQ",std+i).replace("NAME",i)
        print (x,file=sys.stderr)
        os.system(x)
