#!/usr/bin/env python3

import sys
import os
import tree_reader
import seq
import math
import subprocess
import argparse as ap

TEMPNAME = "TEMP"
TEMPCOUNT = 0
scale_trees_to_1 = True

class AlnSet:
    def __init__(self):
        self.alns = set()
        self.filename = ""
        self.lastweight = 0
        self.lastfailweight = None
        self.aicc = 0
        self.aic = 0
        self.bic = 0
        self.like = 0
        self.tree = ""
    
    def __str__(self):
        return str(",".join(list(self.alns)))+" mw:"+str(self.lastweight)+" l:"+str(self.like)+" a:"+str(self.aicc)+" b:"+str(self.bic)+" f:"+str(self.filename)

    def str_for_file(self):
        st = "===\n"
        st += "genes: "+",".join(list(self.alns))+"\n"
        st += "likelihood: "+str(self.like)+"\n"
        st += "aicc: "+str(self.aicc)+"\n"
        st += "bic: "+str(self.bic)+"\n"
        st += "filename: "+str(self.filename)+"\n"
        st += "tree: "+str(self.tree)+"\n"
        return st

def get_lk_from_iqtree(fn):
    of = open(fn+".iqtree","r")
    ll, aic, aicc = 0,0,0
    for j in of:
        if "Log-likelihood of the tree" in j:
            sc = float(j.strip().split(": ")[1].split(" ")[0])
            ll = sc
        elif "Akaike information criterion (AIC) score" in j:
            sc = float(j.strip().split(": ")[1].split(" ")[0])
            aic = sc
        elif "Corrected Akaike information criterion (AICc) score" in j:
            sc = float(j.strip().split(": ")[1].split(" ")[0])
            aicc = sc
        elif "Bayesian information criterion (BIC) score" in j:
            sc = float(j.strip().split(": ")[1].split(" ")[0])
            bic = sc
    of.close()
    return ll, aic, aicc, bic

def concat_nms(alnset1,alnset2,di):
    global TEMPCOUNT
    n = " ".join([di+i for i in list(alnset1.alns)])
    n += " "+" ".join([di+i for i in list(alnset2.alns)])
    tfn = TEMPNAME+str(TEMPCOUNT)
    cmd = "pxcat -s "+n+" -o "+di+tfn+" -p "+di+tfn+".parts"
    os.system(cmd)
    TEMPCOUNT += 1
    print ("  "+cmd,file=sys.stderr)
    return tfn

def run_iqtree(fn,di,nthreads):
    cmd = "iqtree -nt "+nthreads+" -s "+di+fn+" -m GTR+G -pre "+di+fn+" -redo -bb 1000 >> iqtreelog"
    print ("  "+cmd,file=sys.stderr)
    os.system(cmd)
    os.remove(di+fn+".ckp.gz")
    os.remove(di+fn+".log")
    os.remove(di+fn+".mldist")
    os.remove(di+fn+".bionj")

def run_iqtree_part(fn,di,edges,nthreads):
    cmd = "iqtree -nt "+nthreads+" -s "+di+fn+" -m GTR+G -pre "+di+fn+" -"+edges+" "+di+fn+".parts -bb 1000 -redo >> iqtreelog"
    print ("  "+cmd,file=sys.stderr)
    os.system(cmd)
    os.remove(di+fn+".ckp.gz")
    os.remove(di+fn+".log")
    os.remove(di+fn+".mldist")
    os.remove(di+fn+".bionj")

def read_rf_read_to_graph(filename,trmap,di,threads):
    rf = open(filename,"r")
    import networkx as nx
    import numpy as np
    G = nx.Graph()
    sortededges = {}#key is weight, list 
    weights = []
    sets = {} #key is aln , value is the AlnSet
    checked = {} # key is alnset, value is list of alnset
    for i in rf:
        spls = i.strip().split(":")
        t = spls[0].split()
        f1 = trmap[t[0]]
        f2 = trmap[t[1]]
        if f1 not in sets:
            x = AlnSet()
            x.alns.add(f1)
            fn = f1
            # gene trees need to be run to make the trees so no need to redo it
            # run_iqtree(fn,di)
            num,aic,aicc,bic = get_lk_from_iqtree(di+fn)
            x.like = num
            x.aicc = aicc
            x.aic = aic
            x.bic = bic
            x.filename = fn
            x.tree = open(di+fn+".treefile","r").readline()
            sets[f1] = x
            checked[x] = set()
        if f2 not in sets:
            x = AlnSet()
            x.alns.add(f2)
            fn = f2
            # gene trees need to be run to make the trees so no need to redo it
            # run_iqtree(fn,di)
            num,aic,aicc,bic = get_lk_from_iqtree(di+fn)
            x.like = num
            x.aicc = aicc
            x.aic = aic
            x.bic = bic
            x.filename = fn
            x.tree = open(di+fn+".treefile","r").readline()
            sets[f2] = x
            checked[x] = set()
        t = spls[1].strip().split()
        v = float(t[0])
        try:
            sortededges[v].append((f1,f2))
        except:
            sortededges[v] = []
            sortededges[v].append((f1,f2))
            weights.append(v)
        G.add_edge(f1,f2,weight=v)
    rf.close()
    weights.sort()
    return G,weights,sortededges,sets,checked

def run_clustering(weights,sortededges,sets,checked,di,usebic,useaic,edges,stop,stopweight,nthreads):
    going = True
    tempfiles = []
    for i in weights:
        print ("weight:",i,file=sys.stderr)
        if stop: # stop early
            if i >= stopweight:
                print ("STOPPING BECAUSE STOP WAS ENTERED FOR",stopweight,file=sys.stderr)
                break
        for j in sortededges[i]:
            #print >> sys.stderr, sets[j[0]],sets[j[1]]
            if sets[j[0]] == sets[j[1]]:
                continue
            if sets[j[1]] in checked[sets[j[0]]]:
                continue
            if sets[j[0]] in checked[sets[j[1]]]:
                continue
            print ("  ",j,file=sys.stderr)
            if len(sets[j[0]].alns.intersection(sets[j[1]].alns)) > 0:
                print ("SOMETHING WRONG WITH OVERLAP",file=sys.stderr)
                sys.exit(0)
            lks = sets[j[0]].like + sets[j[1]].like
            if usebic:
                ais = sets[j[0]].bic + sets[j[1]].bic
            elif useaic:
                ais = sets[j[0]].aic + sets[j[1]].aic
            else:
                ais = sets[j[0]].aicc + sets[j[1]].aicc
            tfn = concat_nms(sets[j[0]],sets[j[1]],di)
            tempfiles.append(tfn)
            fail = False
            run_iqtree_part(tfn,di,edges,nthreads)
            ll,aic,aicc,bic = get_lk_from_iqtree(di+tfn)
            if usebic:
                comp = bic
            elif useaic:
                comp = aic
            else:
                comp = aicc
            if comp - ais < 0:
                print ("oh yeah!",file=sys.stderr)
                alset = set()
                sets[j[0]].alns = sets[j[0]].alns.union(sets[j[1]].alns)
                for k in sets[j[0]].alns:
                    sets[k] = sets[j[0]]
                for k in sets[j[1]].alns:
                    sets[k] = sets[j[0]]
                #sets[j[1]] = sets[j[0]]
                sets[j[0]].lastweight = i
                sets[j[0]].like = ll
                sets[j[0]].aicc = aicc
                sets[j[0]].aic = aic
                sets[j[0]].bic = bic
                sets[j[0]].filename = tfn
                sets[j[0]].tree = open(di+tfn+".treefile","r").readline()
                #print "++",sets[j[0]]
            else:
                fail = True
                sets[j[0]].lastfailweight = i
                sets[j[1]].lastfailweight = i
                #delete files
                os.remove(di+tfn)
                os.remove(di+tfn+".parts")
                os.remove(di+tfn+".treefile")
                os.remove(di+tfn+".iqtree")
            checked[sets[j[0]]].add(sets[j[1]])
            checked[sets[j[1]]].add(sets[j[0]])
            #x = sorted(G[sortededges[i][0][0]].items(), key=lambda edge: edge[1]['weight'])
    return tempfiles

def make_trees(di,filenamepat,nthreads):
    fls = []
    count = 0
    trmap = {}
    trmapfile = open("treemap","w")
    for i in os.listdir(di):
        if filenamepat in i:
        #if i.split(".")[-1] in seqfileendings:
            run_iqtree(i,di,nthreads)
            fls.append(di+i+".treefile")
            trmap[str(count)] = i
            trmapfile.write(i+"\n")
            count += 1
    trmapfile.close()
    mltreefile = di+"mltrees"
    cmd = "cat "+" ".join(fls)+" >"+ mltreefile
    os.system(cmd)
    return mltreefile,trmap

def scale_treelength(tree,length):
    s = 0
    for j in tree.iternodes():
        s += j.length
    v = (1./s)
    for j in tree.iternodes():
        j.length = j.length * v
    return tree

def calc_rfw(treefile):
    rffile = "rfw"
    if scale_trees_to_1:
        ntreefile = treefile+"_sc"
        ntf = open(ntreefile,"w")
        for i in tree_reader.read_tree_file_iter(treefile):
            ntf.write(scale_treelength(i,1.0).get_newick_repr(True)+";\n")
        ntf.close()
        treefile = ntreefile
    cmd = "bp -t "+treefile+" -rfwp -v > "+rffile+" 2> bplog"
    os.system(cmd)
    return rffile

def read_treemap(filename):
    trmap = {}
    tm = open(filename,"r")
    count = 0 
    for i in tm:
        trmap[str(count)] = i.strip()
        count += 1
    tm.close()
    return trmap

def generate_argparser():
    parser = ap.ArgumentParser(prog="test_clusters.py",
        formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d","--dir",type=str,required=False,default="./",
        help=("Directory with the genes (and trees)."))
    parser.add_argument("-m","--treemap",type=str,required=False,
        help=("List of ordered tree files (each line would be a gene name - for rfw)"))
    parser.add_argument("-w","--rfw",type=str,required=False,
        help=("The RFW file with numbers corresponding to the lines in treemap file."))
    parser.add_argument("-e","--edge",type=str,required=True,
        help=("Edge type (q,spp, or sp)"))
    parser.add_argument("-b","--bic",action='store_true',required=False,default=False,
        help=("Use BIC instead of AICc."))
    parser.add_argument("-a","--aic",action='store_true',required=False,default=False,
        help=("Use AIC instead of AICc."))
    parser.add_argument("-s","--stop",type=float,required=False,
        help=("Stop weight."))
    parser.add_argument("-t","--threads",type=str,required=False,default="2",
        help=("Number of threads."))
    parser.add_argument("-p", "--filenamepattern",help="is there a particular common part of the name (e.g., 'fa')?",
                        default="",required =False )
    return parser

def main():
    parser = generate_argparser()
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")
    args = parser.parse_args(sys.argv[1:])
    di = args.dir
    if di[-1] != "/":
        di += "/"
    rffile = None
    if args.treemap and args.rfw:
        print  ("reading treemap and rfw",file=sys.stderr)
        trmap = read_treemap(args.treemap)
        rffile = args.rfw
    else: # calculate the trees
        print ("making ml trees for genes",file=sys.stderr)
        mltreesfile,trmap = make_trees(di)
        print ("calculating rfw",file=sys.stderr)
        rffile = calc_rfw(mltreesfile)
    
    stop = False
    stopweight = 0
    if args.stop != None:
        stopweight = args.stop
        stop = True

    print ("creating graph",file=sys.stderr)
    G,weights,sortededges,sets,checked = read_rf_read_to_graph(rffile,trmap,di,args.threads)
    
    usebic = False
    useaic = False
    if args.bic == True:
        print ("using BIC",file=sys.stderr)
        usebic = True
    elif args.aic == True:
        print ("using AIC",file=sys.stderr)
        useaic = True
    edges = args.edge
    
    print ("clustering",file=sys.stderr)
    tempfiles = run_clustering(weights,sortededges,sets,checked,di,usebic,useaic,edges,stop,stopweight,args.threads)

    print ("writing final sets to finalsets.txt",file=sys.stderr)
    outfile = open("finalsets.txt","w")
    print ("")
    print ("FINAL SETS")
    print ("==========")
    keepers = set()
    for i in list(set(sets.values())):
        print (i)
        outfile.write(i.str_for_file()+"\n")
        keepers.add(i.filename)
    outfile.close()

    print ("deleting other files")
    for i in tempfiles:
        if i not in keepers:
            if os.path.isfile(di+i):
                os.remove(di+i)
                os.remove(di+i+".parts")
                os.remove(di+i+".treefile")
                os.remove(di+i+".iqtree")
            

if __name__ == "__main__":
    main()
