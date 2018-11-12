import sys
import os
import os.path

cmd = "iqtree -nt 3 -s SEQ -alrt 1000 -sp SEQ.parts -pre NAME -m GTR+G -redo"
cmd = "iqtree -nt 3 -s SEQ -alrt 1000 -pre NAME -m GTR+G -redo"

if __name__ == "__main__":
    if len (sys.argv) != 2 and len(sys.argv) != 3:
        print "python "+sys.argv[0]+" dir (inseqname)" 
        sys.exit(0)
    st = sys.argv[1]
    isn = ""
    if len(sys.argv) == 3:
        isn = sys.argv[2]
    if st[-1] != "/":
        st += "/"
    for i in os.listdir(st):
        if len(isn) > 0:
            if isn not in i:
                continue
        x = cmd.replace("SEQ",st+i).replace("NAME",i)
        print x
        os.system(x)
