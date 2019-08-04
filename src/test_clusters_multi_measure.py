#!/usr/bin/env python3

import sys
import os
import tree_reader
import seq
import math
import subprocess
import argparse as ap

"""
this is very similar to the test_clusters.py one just that it only compares things if they 
meet certain criteria like RFP <= some value
"""

TEMPNAME = "TEMP"
TEMPCOUNT = 0
scale_trees_to_1 = True

# k = num params
# n = num of sites
def calc_aicc(aic,n,k):
    x = aic + (((2*math.pow(k,2))+(2*k))/(n-k-1))
    return x

from test_clusters import get_lk_from_iqtree,concat_nms,run_iqtree,run_iqtree_part
from test_clusters import make_trees,scale_treelength,calc_rfw,read_treemap
from test_clusters import AlnSet


def read_rf_read_to_graph(filenamerfwpv,filenamerfp,trmap,di):
    rf = open(filenamerfwpv,"r")
    rfp = open(filenamerfp,"r")
    rfpmap = {} # 
    for i in rfp:
        spls = i.strip().split(":")
        t = spls[0].split()
        f1 = trmap[t[0]]
        f2 = trmap[t[1]]
        if f1 not in rfpmap:
            rfpmap[f1] = {}
        if f2 not in rfpmap:
            rfpmap[f2] = {}
        rfpmap[f1][f2] = float(spls[1].strip())
        rfpmap[f2][f1] = float(spls[1].strip())
    rfp.close()
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
        m = float(t[1])
        try:
            sortededges[v].append((f1,f2))
        except:
            sortededges[v] = []
            sortededges[v].append((f1,f2))
            weights.append(v)
        G.add_edge(f1,f2,rfwp=v,maxedgediff=m,rfp=rfpmap[f1][f2])
    rf.close()
    weights.sort()
    return G,weights,sortededges,sets,checked

def run_clustering(weights,sortededges,sets,checked,di,usebic,useaic,edges,stop,stopweight,G):
    going = True
    tempfiles = []
    for i in weights:
        print ("weight:",i,file=sys.stderr)
        if stop: # stop early
            if i >= stopweight:
                print ("STOPPING BECAUSE STOP WAS ENTERED FOR",stopweight,file=sys.stderr)
                break
        for j in sortededges[i]:
            if G[j[0]][j[1]]['rfp'] > 0:
                continue
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

def generate_argparser():
    parser = ap.ArgumentParser(prog="test_clusters_multi_measure.py",
        formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d","--dir",type=str,required=False,default="./",
        help=("Directory with the genes (and trees)."))
    parser.add_argument("-m","--treemap",type=str,required=True,
        help=("List of ordered tree files (each line would be a gene name - for rfw)"))
    parser.add_argument("-w","--rfwpv",type=str,required=True,
        help=("The RFWPV file with numbers corresponding to the lines in treemap file. bp -t t -rfwp -v."))
    parser.add_argument("-r","--rfp",type=str,required=True,
        help=("The RFP file with numbers corresponding to the lines in treemap file. bp -t t -rfp -scut 90 -v."))
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
    print ("reading treemap and rfw",file=sys.stderr)
    trmap = read_treemap(args.treemap)
    rffile = args.rfwpv
    rfpfile = args.rfp
    
    stop = False
    stopweight = 0
    if args.stop != None:
        stopweight = args.stop
        stop = True

    print ("creating graph",file=sys.stderr)
    G,weights,sortededges,sets,checked = read_rf_read_to_graph(rffile,rfpfile,trmap,di,args.threads)
    
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
    tempfiles = run_clustering(weights,sortededges,sets,checked,di,usebic,useaic,edges,stop,stopweight,G,args.threads)

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
