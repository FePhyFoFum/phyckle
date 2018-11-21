import sys
import os
"""
this is meant to take the constraint file
tree
-tree
-tree

and output the csv files with whether the mltrees match the constraints
"""

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "python "+sys.argv[0]+" constraints mltreedir"
        sys.exit(0)
    mldir = sys.argv[2]
    tf = open(sys.argv[1],"r")
    count = 0
    for i in tf:
        if i[0] == "-":
            i = i[1:].strip()
        else:
            i = i.strip()
        of = open("constraints"+str(count),"w")
        of.write(i+';')
        of.close()
        cmd = "python ~/Dropbox/programming/phyckle/src/check_outgroup_mono.py constraints"+str(count)+" "+mldir+" > constraints"+str(count)+".csv"
        print cmd
        os.system(cmd)
        count += 1
    tf.close()
