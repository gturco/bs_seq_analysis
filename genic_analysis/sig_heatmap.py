from collections import defaultdict
import scipy.stats.mstats as s
from numpy import array

def main(mbin,meth_data):
    d = meth_dict(meth_data)
    mmax,mtypes, mtype = stats_points(d)
    scores(mmax,mbin,mtypes,mtype,d)



def meth_dict(meth_data):
    d = defaultdict(list)
    for line in open(meth_data):
        pos,fract,mtype = line.strip().split("\t")
        d[mtype].append((int(pos),float(fract)))
    return d

def stats_points(m):
    mtypes = []

    for mtype in m.keys():
        ###mtypes.append(sorted(d[mtype], key=lambda x:x[0]))
        meth = [f for p,f in sorted(m[mtype]) ]
        mtypes.append(meth)

    mmax = len(mtypes[0])
    return mmax, mtypes, mtype

def scores(mmax,mbin,mtypes,mtype,m):
    for n in range(0,mmax):
        a = mtypes[0][n:n+mbin]
        b = array(mtypes[1][n:n+mbin])
        c = array(mtypes[2][n:n+mbin])
        d = array(mtypes[3][n:n+mbin])
        hstat, pval = s.kruskalwallis(a,b,c,d)
        if n < 100:
            print "{0}\t{1}".format(sorted(m[mtype])[n][0] ,pval)
            print "{0}\t{1}".format(sorted(m[mtype])[n][0] + 10 ,pval)
        elif n > 199:
            print "{0}\t{1}".format(sorted(m[mtype])[n][0] - 100,pval)
            print "{0}\t{1}".format(sorted(m[mtype])[n][0] - 90 ,pval)
        else:
            print "{0}\t{1}".format(sorted(m[mtype])[n][0] ,pval)


### 10,20,30,-1000 (100 points)
    ### 20,40 -10 points

main(20,"/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CG_all_4.txt")


###ggplot(x, aes(pos, "test")) + geom_tile(aes(fill = sig))

