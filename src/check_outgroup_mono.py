import sys
import os
import subprocess

"""
this will use bp to check for whether the outgroup is monophyletic

it makes some assumptions about files. the ml files are just X.fasta.treefile
"""
#front = "RAxML_bipartitions"
tail = ".treefile"

# SH, 
BPCMD = "bp -t MLTREEFILE -c OUTGROUP -scut 80"

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "python "+sys.argv[0]+" outgroup.tre mldir"
        sys.exit(0)
    
    outgroupf = sys.argv[1]
    print "gene,"+outgroupf
    gened = sys.argv[2]
    if gened[-1] != "/":
        gened += "/"
    for i in os.listdir(gened):
        if tail in i:
        #if front in i:
            gf = i
            cmd = BPCMD.replace("MLTREEFILE",gened+i).replace("OUTGROUP",outgroupf)
            #print cmd
            process = subprocess.Popen(cmd.split(" "),stdout = subprocess.PIPE,stderr=subprocess.PIPE)
            out,err = process.communicate()
            count = 0
            start = False
            for j in out.split("\n"):
                if len(j.strip()) > 10:
                    count += 1
            print gf+","+str(count)
