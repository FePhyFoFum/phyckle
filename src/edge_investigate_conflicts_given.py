#!/usr/bin/env python3

import argparse
import tree_reader
import sys
import os
import seq

# do we want to have a minimum branch length below which you ignore?
MINBL = 0.0000

RAXMLLINE = "RAXML -T THREAD -f a -# 100 -x 1234 -s SEQ -n NAME -p 1234 -m GTRCAT -g CONFILE"
RAXMLLINE = "RAXML -T THREAD -s SEQ -n NAME -p 1234 -m GTRCAT -g CONFILE"
run = True # can be helpful if you want to distribute runs

class Bipart:
    def __init__ (self,lf,rt):
        self.left = lf
        self.right = rt
        self.union = lf.union(rt)

    def __str__(self):
        x = ",".join(list(self.left))
        y = ",".join(list(self.right))
        return x+" | "+y

    def conflict(self, inbp):
        if len(inbp.right.intersection(self.right)) > 0 and len(inbp.right.intersection(self.left)) > 0:
            if len(inbp.left.intersection(self.right)) > 0 and len(inbp.left.intersection(self.left)) > 0 :
                return True
        if len(inbp.left.intersection(self.left)) > 0 and len(inbp.left.intersection(self.right)) > 0:
            if len(inbp.right.intersection(self.left)) > 0 and len(inbp.right.intersection(self.right)) > 0:
                return True
        return False
    
    def equal(self, inbp):
        inter = self.union.intersection(inbp.union)
        if self.left.intersection(inter) == inbp.left.intersection(inter)  and self.right.intersection(inter)  == inbp.right.intersection(inter) :
            return True
        if self.left.intersection(inter)  == inbp.right.intersection(inter) and self.right.intersection(inter)  == inbp.left.intersection(inter) :
            return True
        return False

    def newick(self, exclude=None):
        if exclude == None:
            return "(("+",".join(list(self.left))+"),"+",".join(list(self.right))+")"
        else:
            if len(self.left.difference(exclude)) < 2:
                return "DONTRUN"
            elif len(self.right.difference(exclude)) < 2:
                return "DONTRUN"
            else:
                return "(("+",".join(list(self.left.difference(exclude)))+"),"+",".join(list(self.right.difference(exclude)))+")"

def get_biparts(rt):
    bps = []
    out = set(rt.lvsnms())
    for i in rt.iternodes():
        if len(i.children) == 0:
            continue
        if i == rt:
            continue
        if MINBL > 0:
            if rt.length < MINBL:
                continue
        right = set(i.lvsnms())
        left = set(out-right)
        bp = Bipart(right,left)
        bps.append(bp)
    return bps

qsubstring = """#PBS -A eebsmith_flux
#PBS -q flux
#PBS -M eebsmith@umich.edu
#PBS -m MAILOPTION
#PBS -j oe
#PBS -V
#PBS -N raxml.RUN
#PBS -l nodes=1:ppn=8,mem=4gb,walltime=168:00:00
cd "$PBS_O_WORKDIR"
COMMAND
"""
def create_qsub(count,outfilename,runname,command):
    # MAILOPTION should be a abort, b begin, e end, you can do abe, ab, ae, be, a, b, e
    if count % 10 == 0:
        mo = "ae"
    else:
        mo = "a"
    of = open(outfilename,"w")
    of.write(qsubstring.replace("COMMAND",command).replace("MAILOPTION",mo).replace("RUN",runname))
    of.close()


"""
conflicting bipartitions should have a - in the front

"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-c","--constraint",help="What is the file with the constraints",required=True)
    parser.add_argument("-o","--outfile",help="What is the outfile",required=True)
    parser.add_argument("-m","--outdir",help="What is the outdir for all the files",required=True)
    parser.add_argument("-r", "--raxml", help="Where is the raxml bin?",default="raxml",required=False)
    parser.add_argument("-t", "--threads",help="How many threads for raxml?",default="4",required=False)
    parser.add_argument("-d", "--indir",help="Which directory?",required=True)
    parser.add_argument("-s", "--filenamepattern",help="is there a particular common part of the name?",
                        default="",required =False )

    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    args = parser.parse_args()

    isn = args.filenamepattern

    RAXMLLINE = RAXMLLINE.replace("THREAD",args.threads).replace("RAXML",args.raxml)

    test_bps = {}
    test_bp_con_sets = {} #key is int of testbipart and value is list of unique conflicts
    outf = open(args.outfile,"w")
    count = 0
    curcount = 0
    of = open(args.constraint,"r")
    for i in of:
        tree = None
        if i[0] == "-":
            tree = tree_reader.read_tree_string(i[1:])
            if curcount not in test_bp_con_sets:
                test_bp_con_sets[curcount] = []
            test_bp_con_sets[curcount].append(get_biparts(tree)[0]) # only adding the first one
        else:
            tree = tree_reader.read_tree_string(i)
            test_bps[count] = get_biparts(tree)[0] # only getting the first bp
            curcount = count
            count += 1
    of.close()
    
    # write the constraints to a file
    for i in test_bp_con_sets:
        outf.write("constraint: "+str(i)+"\n")
        outf.write(str(test_bps[i])+"\n")
        print ("constraint:",i,file=sys.stderr)
        print (" -conflicts-",file=sys.stderr)
        outf.write(" -conflicts-\n")
        for j in range(len(test_bp_con_sets[i])):
            print (" ",j,file=sys.stderr)
            outf.write(" "+str(j)+" : "+str(test_bp_con_sets[i][j])+"\n")
    outf.close()
    #sys.exit(0)
    # run all the raxml things
    print ("running constraints",file=sys.stderr)
    gdir = args.indir
    if gdir[-1] != "/":
        gdir += "/"
    odir = args.outdir
    if os.path.isdir(odir) == False:
        print ("making",odir,file=sys.stderr)
        os.mkdir(odir)
    if odir[-1] != "/":
        odir += "/"
    count = 0
    for i in os.listdir(gdir):
        seqf = i
        print ("SEQUENCE:",seqf,file=sys.stderr)
        seq_names = set()
        for j in seq.read_fasta_file_iter(gdir+i):
            seq_names.add(j.name)
        for j in test_bp_con_sets:
            print (" constraint:",j,file=sys.stderr)
            ex = seq_names.symmetric_difference(test_bps[j].union)
            constring = ""
            if len(ex) > 0:
                constring = test_bps[j].newick(ex)
                if constring == "DONTRUN":
                    print ("DONTRUN",file=sys.stderr)
                    continue
            else:
                constring = test_bps[j].newick()
            confilename = odir+i+"___cons_"+str(j)
            cf = open(confilename,"w")
            cf.write(constring+";")
            cf.close()
            print (" ",constring,file=sys.stderr)
            name = i+"___cons_"+str(j)
            cmd = RAXMLLINE.replace("SEQ",gdir+i).replace("NAME",name).replace("CONFILE",confilename)
            print (cmd,file=sys.stdout)
            if run == True:
                os.system(cmd)
                os.system("mv RAxML_bipartitions."+name+" "+odir)
                os.system("mv RAxML_info."+name+" "+odir)
                os.system("rm RAxML_bipartitionsBranchLabels."+name)
                os.system("rm RAxML_bootstrap."+name)
                os.system("rm RAxML_bestTree."+name)
            print ("  -conflicts-",file=sys.stderr)
            for k in range(len(test_bp_con_sets[j])):
                ex = seq_names.symmetric_difference(test_bp_con_sets[j][k].union)
                print (ex,seq_names,test_bp_con_sets[j][k].union,file=sys.stderr)
                constring = ""
                if len(ex) > 0:
                    constring = test_bp_con_sets[j][k].newick(ex)
                    if constring == "DONTRUN":
                        print ("DONTRUN",file=sys.stderr)
                        continue
                else:
                    constring = test_bp_con_sets[j][k].newick()
                print ("  ",k,file=sys.stderr)
                confilename = odir+i+"___cons_"+str(j)+"_conf_"+str(k)
                cf = open(confilename,"w")
                cf.write(constring+";")
                cf.close()
                print ("  ",constring,file=sys.stderr)
                name = i+"___cons_"+str(j)+"_conf_"+str(k)
                cmd = RAXMLLINE.replace("SEQ",gdir+i).replace("NAME",name).replace("CONFILE",confilename)
                print (cmd,file=sys.stdout)
                if run == True:
                    os.system(cmd)
                    os.system("mv RAxML_bipartitions."+name+" "+odir)
                    os.system("mv RAxML_info."+name+" "+odir)
                    os.system("rm RAxML_bipartitionsBranchLabels."+name)
                    os.system("rm RAxML_bootstrap."+name)
                    os.system("rm RAxML_bestTree."+name)
