#!/usr/bin/env python3

import sys
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--csvfiles", nargs='+',help="List each file in the order to process.",required=True)
    parser.add_argument("-o","--outfile",help="Outfile for analysis",required=True)
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")
    args = parser.parse_args()
    print(args.csvfiles)

    nms = {}
    files = args.csvfiles
    fileheaders = []
    count = 0
    outf=open(args.outfile,"w")
    for i in files:
        found = set()
        fl = open(i,"r")
        header = fl.readline().strip().split(",")
        fileheaders.append(header[1:])
        add = 0
        for j in fl:
            j = j.strip().split(",")
            # can edit the first thing in case there are directories and such, but probably better to sed those
            try:
                nms[j[0]]
            except:
                nms[j[0]] = []
                for k in range(count):
                    nms[j[0]].append("-")
            found.add(j[0])
            for k in j[1:]:
                nms[j[0]].append(k)
            add = len(j[1:])
        fl.close()
        notfound = set(nms.keys())-found
        for j in notfound:
            for k in range(add):
                nms[j].append("-")
        count += add
    stri = "gene"
    for i in fileheaders:
        stri += ","+",".join(i)
    print (stri,file=outf)
    for i in nms:
        print (i+","+",".join(nms[i]),file=outf)
    outf.close()
