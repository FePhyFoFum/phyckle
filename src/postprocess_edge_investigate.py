import sys
import os

"""
out.conflict
"""

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "python "+sys.argv[0]+" out.conflict outd"
        sys.exit(0)
    oc = open(sys.argv[1],"r")
    constraints = {} #key is constraint number and value is conflict alts
    curcons = None
    curlist = None
    reccon = False
    for i in oc:
        if "constraint" in i:
            reccon = False
            spls = i.strip().split(": ")[1]
            if curcons == None:
                curcons = spls
                curlist = []
            else:
                constraints[curcons] = curlist
                curlist = []
                curcons = spls
            continue
        if "-conflicts-" in i:
            reccon = True
            continue
        if reccon == True:
            spls = i.strip().split(" : ")[0]
            curlist.append(spls)
    #get last one
    constraints[curcons] = curlist
    oc.close()
    outd = sys.argv[2]
    if outd[-1] != "/":
        outd += "/"
    
    scores = {} #key is seqfname , value is dictionary with runname
    analyses = set ()
    for i in os.listdir(outd):
        if "info" in i:
            spl = i.split("RAxML_info.")[-1].split("___")
            sfn = spl[0]
            if sfn not in scores:
                scores[sfn] = {}
            anafn = spl[1]
            analyses.add(anafn)
            # print sfn,anafn
            of = open(outd+i,"r")
            for j in of:
                if "Final GAMMA-based" in j:
                    sc = float(j.strip().split()[-1])
                    scores[sfn][anafn] = sc
            of.close()

    
    for i in constraints:
        outf = open("cons_"+i+".csv","w")
        nm = "cons_"+i
        outf.write("gene,"+nm)
        for j in constraints[i]:
            nm2 = ","+nm+"_conf_"+str(j)
            outf.write(nm2)
        outf.write(",bestone,best,secondbest,diffbestsecondbest\n")
        print nm
        for k in scores:
            stri = k +","+str(scores[k][nm])
            best = scores[k][nm]
            secondbest = -9999999999
            bestone = nm
            for j in constraints[i]:
                nm2 = nm+"_conf_"+str(j)
                try:
                    stri += ","+str(scores[k][nm2])
                    if scores[k][nm2] > best:
                        secondbest = best
                        best = scores[k][nm2]
                        bestone = nm2
                        continue
                    if secondbest == -9999999999:
                        secondbest = scores[k][nm2]
                    else:
                        if scores[k][nm2] > secondbest and scores[k][nm2] < best:
                            secondbest = scores[k][nm2]
                except:
                    stri += ",-"

            stri += ","+bestone+","+str(best)+","+str(secondbest)+","+str(best-secondbest)
            outf.write(stri+"\n")
        outf.close()
            
        #for j in scores:
        #    stri = j
    """
    for i in constraints:
        outf = open("cons_"+i+".out","w")
        nm = "cons_"+i
        print nm
        for j in scores:
            stri = str(j)+","+str(scores[j][nm])
            #for k in constraints[i]:
                #print j, scores[i][j]
            #    print j
        outf.close()
    """
