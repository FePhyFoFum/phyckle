import sys

"""
this is intended to process the combined files from doing a combine_csv
there is some variation with these files but in general, it expects there to be a 
constraints#
and then a bunch of things from the basic cons_0.csv
"""

verbose = True

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "python "+sys.argv[0]+" combined.csv"
        sys.exit(0)
    fl = open(sys.argv[1],"r")
    spls = fl.readline().strip().split(",") # read the first line
    genei = 0
    bestone = None
    diff = None
    constraintscount = 0
    constraintd = {}
    constraintdifs = {}
    constraintgenes = {}
    constraint2difs = {}
    constraint2genes = {}
    constraintdifsbad = {}
    constraintgenesbad = {}
    constraint2difsbad = {}
    constraint2genesbad = {}
    requiredi = None
    translation = {}
    badoutgroupi = None
    badblleni = None
    badoutgroup = set() # these are the ones that are bad from outgroup if there is an outgroup in the list of things
    badbllen = set() # these are the ones that are bad from bl if there is a bad max bl
    count = 0
    for i in spls:
        if i == "bestone":
            bestone = count
        if i == "diffbestsecondbest":
            diff = count
        if "constraint" in i:
            constraintd[constraintscount] = count
            constraintgenes[constraintscount] = []
            constraintdifs[constraintscount] = []
            constraint2genes[constraintscount] = []
            constraint2difs[constraintscount] = []
            #bad ones
            constraintgenesbad[constraintscount] = []
            constraintdifsbad[constraintscount] = []
            constraint2genesbad[constraintscount] = []
            constraint2difsbad[constraintscount] = []
            if constraintscount == 0:
                translation["cons_0"] = constraintscount
            else:
                translation["cons_0_conf_"+str(constraintscount-1)] = constraintscount
            constraintscount += 1
        if "outgroup" in i:
            badoutgroupi = count
        if "maxbl" in i:
            badblleni = count
        if "required" == i:
            requiredi = count
        count += 1
    for i in fl:
        spls = i.strip().split(",")
        if "-" in spls:
            continue
        bad = False
        if badoutgroupi != None:
            if int(spls[badoutgroupi]) > 0:
                badoutgroup.add(spls[0])
                bad = True
        if badblleni != None:
            if float(spls[badblleni]) > 2.5:
                badbllen.add(spls[0])
                bad = True
        if requiredi != None:
            if int(spls[requiredi]) == 1:
                continue
        t = translation[spls[bestone]]
        mlcon = int(spls[constraintd[t]])
        if mlcon == 0:
            #print spls[genei],spls[bestone],mlcon,spls[diff],t
            constraintdifs[t].append(float(spls[diff]))
            constraintgenes[t].append(spls[genei])
            if bad == False:
                constraintgenesbad[t].append(spls[genei])
                constraintdifsbad[t].append(float(spls[diff]))
            if float(spls[diff]) >= 2:
                constraint2difs[t].append(float(spls[diff]))
                constraint2genes[t].append(spls[genei])
                if bad == False:
                    constraint2difsbad[t].append(float(spls[diff]))
                    constraint2genesbad[t].append(spls[genei])
    fl.close()
    print "constraint sumdiff sum2diff lengenes len2genes"
    for i in constraintdifs:
        print i,sum(constraintdifs[i]),sum(constraint2difs[i]),len(constraintgenes[i]),len(constraint2genes[i])
        import numpy
        x = numpy.argsort(constraintdifs[i])[::-1][0:3]
        for j in x:
            print "  ",constraintdifs[i][j],constraintgenes[i][j]


    lks = []
    lks2 = []
    gns = []
    gns2 = []
    totalgn = 0
    totalgn2 = 0
    print "\noutgroup and ml"
    print "constraint sumdiff sum2diff lengenes len2genes"
    for i in constraintdifs:
        print i,sum(constraintdifsbad[i]),sum(constraint2difsbad[i]),len(constraintgenesbad[i]),len(constraint2genesbad[i])
        lks.append(str(sum(constraintdifsbad[i])))
        lks2.append(str(sum(constraint2difsbad[i])))
        gns.append(str(len(constraintgenesbad[i])))
        gns2.append(str(len(constraint2genesbad[i])))
        totalgn += len(constraintgenesbad[i])
        totalgn2 += len(constraint2genesbad[i])
        import numpy
        x = numpy.argsort(constraintdifsbad[i])[::-1][0:3]
        for j in x:
            print "  ",constraintdifsbad[i][j],constraintgenesbad[i][j]
        #print "cp "+" ".join(constraint2genesbad[i])
    gns.append(str(852-totalgn))
    gns2.append(str(852-totalgn2))

    # write the pies file
    of = open("pies.R","w")
    of.write('svg("pies.svg")\n')
    of.write('library(colortools)\n')
    of.write('par(mfrow=c(2,2))\n')
    of.write('a = c('+",".join(lks)+')\n')
    of.write('b = c('+",".join(lks2)+')\n')
    of.write('d = c('+",".join(gns)+')\n')
    of.write('e = c('+",".join(gns2)+')\n')
    of.write('pie(a,col=splitComp("steelblue",plot=F),labels="")\n')
    of.write('pie(b,col=splitComp("steelblue",plot=F),labels="")\n')
    of.write('pie(d,col=splitComp("steelblue",plot=F),labels="")\n')
    of.write('pie(e,col=splitComp("steelblue",plot=F),labels="")\n')
    of.write('dev.off()')
    of.close()
