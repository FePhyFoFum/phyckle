import sys


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "print "+sys.argv[0]+" files.csv"
        sys.exit(0)
    nms = {}
    files = sys.argv[1:]
    fileheaders = []
    count = 0
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
    print stri
    for i in nms:
        print i+","+",".join(nms[i])