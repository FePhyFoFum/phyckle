#!/usr/bin/env python3

import sys
import os
import os.path
import argparse
import re

"""

probably
rm RAxML_info* RAxML_log* RAxML_fastTree.* RAxML_parsimonyTree* RAxML_result*

"""

cmd = "RAXML -T THREAD -s SEQ -n NAME -p 1234 -m GTRCAT"
cmdsh = "RAXML -T THREAD -f J -s SEQ -n NAME -p 1234 -m GTRCAT -t TREEFILE"

def fix_raxml(infile):
    p = re.compile(r'(:[0-9]\.[0-9]+)\[([0-9]+)\]')
    inf = open(infile,"r")
    strg = inf.readline().strip()
    inf.close()
    result = p.sub(r'\2\1',strg)
    ouf = open(infile,"w") 
    ouf.write(result+"\n")
    ouf.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-r", "--raxml", help="Where is the raxml bin?",default="raxml",required=False)
    parser.add_argument("-t", "--threads",help="How many threads for raxml?",default="4",required=False)
    parser.add_argument("-d", "--dir",help="Which directory?",required=True)
    parser.add_argument("-s", "--filenamepattern",help="is there a particular common part of the name?",
                        default="",required =False )
    parser.add_argument("-p","--run_sh",help="Run SH like test with raxml as well",action='store_true',required=False)

    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    args = parser.parse_args()

    std = args.dir
    isn = args.filenamepattern

    if std[-1] != "/":
        std += "/"
    for i in os.listdir(std):
        if len(isn) > 0:
            if isn not in i:
                continue
        treefile = "RAxML_bestTree."+i
        if os.path.exists(treefile) and args.run_sh:
            x = cmdsh.replace("SEQ",std+i).replace("NAME",i).replace("THREAD",args.threads).replace("RAXML",args.raxml).replace("TREEFILE",treefile)
            print (x,file=sys.stderr)
            os.system(x)
            infile = "RAxML_fastTreeSH_Support."+i
            fix_raxml(infile)
        else:
            x = cmd.replace("SEQ",std+i).replace("NAME",i).replace("THREAD",args.threads).replace("RAXML",args.raxml)
            print (x,file=sys.stderr)
            os.system(x)
            if args.run_sh:
                os.system("rm RAxML_info."+i)
                x = cmdsh.replace("SEQ",std+i).replace("NAME",i).replace("THREAD",args.threads).replace("RAXML",args.raxml).replace("TREEFILE",treefile)
                print (x,file=sys.stderr)
                os.system(x)
                infile = "RAxML_fastTreeSH_Support."+i
                fix_raxml(infile)
