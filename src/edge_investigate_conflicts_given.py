import tree_reader
import sys
import os
import seq

# do we want to have a minimum branch length below which you ignore?
MINBL = 0.0000

# RAXMLLINE = "raxml -T 2 -s SEQ -n NAME -p 1234 -m GTRCAT -g CONFILE"
RAXMLLINE = "raxml -T 8 -s SEQ -n NAME -p 1234 -m GTRCAT -g CONFILE"
IQTREELINE = "iqtree -s SEQ -m GTR+G -g CONFILE -nt 2 -pre NAME"
makeqsub = False
raxml = True

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
    if len(sys.argv) != 5:
        print "python "+sys.argv[0]+ " in.bptree outfile dir_for_genes_fa outdir"
        sys.exit(0)

    test_bps = {}
    test_bp_con_sets = {} #key is int of testbipart and value is list of unique conflicts
    outf = open(sys.argv[2],"w")
    count = 0
    curcount = 0
    of = open(sys.argv[1],"r")
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
        print "constraint:",i
        print " -conflicts-"
        outf.write(" -conflicts-\n")
        for j in range(len(test_bp_con_sets[i])):
            print " ",j
            outf.write(" "+str(j)+" : "+str(test_bp_con_sets[i][j])+"\n")
    outf.close()
    #sys.exit(0)
    # run all the raxml things
    print "running constraints"
    gdir = sys.argv[3]
    if gdir[-1] != "/":
        gdir += "/"
    odir = sys.argv[4]
    if os.path.isdir(odir) == False:
        print "making",odir
        os.mkdir(odir)
    if odir[-1] != "/":
        odir += "/"
    count = 0
    for i in os.listdir(gdir):
        if i[0] == ".":
            continue
        if i[-4:] == ".log":
            continue
        if ".reduced" in  i:
            continue
        if ".parstree" in i:
            continue
        if ".uniqueseq" in i:
            continue
        if ".iqtree" in i:
            continue
        if ".treefile" in i:
            continue
        #if "c1c2" not in i:
        #    continue
        seqf = i
        print "SEQUENCE:",seqf
        seq_names = set()
        for j in seq.read_fasta_file_iter(gdir+i):
            seq_names.add(j.name)
        for j in test_bp_con_sets:
            print " constraint:",j
            ex = seq_names.symmetric_difference(test_bps[j].union)
            constring = ""
            if len(ex) > 0:
                constring = test_bps[j].newick(ex)
                if constring == "DONTRUN":
                    print "DONTRUN"
                    continue
            else:
                constring = test_bps[j].newick()
            confilename = odir+i+"___cons_"+str(j)
            cf = open(confilename,"w")
            cf.write(constring+";")
            cf.close()
            print " ",constring
            name = i+"___cons_"+str(j)
            if raxml:
                cmd = RAXMLLINE.replace("SEQ",gdir+i).replace("NAME",name).replace("CONFILE",confilename)
                print " ",cmd
                if makeqsub:
                    create_qsub(count,name+".qsub",name,cmd)
                else:
                    os.system(cmd)
                    os.system("mv RAxML_bestTree."+name+" "+odir)
                    os.system("mv RAxML_info."+name+" "+odir)
                    os.system("rm RAxML_log."+name)
                    os.system("rm RAxML_result."+name)
            else:
                cmd = IQTREELINE.replace("SEQ",gdir+i).replace("CONFILE",confilename).replace("NAME",name)
                print " ",cmd
                if makeqsub:
                    create_qsub(count,name+".qsub",name,cmd)
                else:
                    os.system(cmd)
                    addname = "___cons_"+str(j)
                    os.system("mv "+gdir+name+".treefile "+odir)
                    os.system("mv "+gdir+name+".log "+odir)
                    os.system("mv "+gdir+name+".iqtree "+odir)
                    os.system("rm "+gdir+name+".parstree ")
                    os.system("rm "+gdir+name+".ckp.gz")
                count += 1
            print "  -conflicts-"
            for k in range(len(test_bp_con_sets[j])):
                ex = seq_names.symmetric_difference(test_bp_con_sets[j][k].union)
                print ex,seq_names,test_bp_con_sets[j][k].union
                constring = ""
                if len(ex) > 0:
                    constring = test_bp_con_sets[j][k].newick(ex)
                    if constring == "DONTRUN":
                        print "DONTRUN"
                        continue
                else:
                    constring = test_bp_con_sets[j][k].newick()
                print "  ",k
                confilename = odir+i+"___cons_"+str(j)+"_conf_"+str(k)
                cf = open(confilename,"w")
                cf.write(constring+";")
                cf.close()
                print "  ",constring
                name = i+"___cons_"+str(j)+"_conf_"+str(k)
                if raxml:
                    cmd = RAXMLLINE.replace("SEQ",gdir+i).replace("NAME",name).replace("CONFILE",confilename)
                    print "  ",cmd
                    if makeqsub:
                        create_qsub(count,name+".qsub",name,cmd)
                    else:
                        os.system(cmd)
                        os.system("mv RAxML_bestTree."+name+" "+odir)
                        os.system("mv RAxML_info."+name+" "+odir)
                        os.system("rm RAxML_log."+name)
                        os.system("rm RAxML_result."+name)
                else:
                    cmd = IQTREELINE.replace("SEQ",gdir+i).replace("CONFILE",confilename).replace("NAME",name)
                    print " ",cmd
                    if makeqsub:
                        create_qsub(count,name+".qsub",name,cmd)
                    else:
                        os.system(cmd)
                        #addname = "___cons_"+str(j)+"_conf_"+str(k)
                        os.system("mv "+gdir+name+".treefile "+odir)
                        os.system("mv "+gdir+name+".log "+odir)
                        os.system("mv "+gdir+name+".iqtree "+odir)
                        os.system("rm "+gdir+name+".parstree ")
                        os.system("rm "+gdir+name+".ckp.gz")
                    count += 1
