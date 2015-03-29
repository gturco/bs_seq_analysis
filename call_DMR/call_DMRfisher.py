###100 bp window 50bp iteration
### less than 20 cyosines in window than dicard
### Calculate for each stage... 
###Kruskal wallis test on each window
## FDR 0.05
### Atleast 2 fold change
### concantate overlaping
from intersection import Intersecter, Feature
from flatfeature import Bed
from collections import defaultdict
import scipy.stats as s


def loadintointersect_meth(meth_file):
    query_list = {}
    seen = set()
    for line in open(meth_file):
        args = line.strip().split("\t")
    ### ADD one because bed adds one too number
        seqid = args[0]
        if seqid not in list(query_list.keys()):
                query_list[seqid] = Intersecter()
        start = args[1]
        if line in seen: continue
        seen.add(line)
        ### remove repeats
        if int(args[4]) > 0:
            args[4] = 1
        if int(args[10]) > 0:
            args[10] = 1

        name = "{0}_{1}".format(args[4],args[10])
        query_list[seqid].add_interval(Feature(int(start), int(start),name=name))
    return query_list





def freq(feature_file,window_size,interval,meth_data):
  features = Bed(feature_file)
  for feature in features:
        region = range(int(feature["start"]),int(feature["end"])+1)
        for window_start in region[::interval]:
            window_end = window_start + window_size
            if window_end > region[-1]:
                matches = meth_data[feature['seqid']].find(window_start, region[-1])
            else:
                matches = meth_data[feature['seqid']].find(window_start,window_end)
            if len(matches) < 15 : continue
            kw(matches,feature['seqid'],window_start,window_end)

def kw(matches,seqid,start,stop):
    pop_a = []
    pop_b = []
    
    for m in matches:
        m_info = m.name
        a,b = map(int,m_info.split("_"))
        pop_a.append(a)
        pop_b.append(b)
        ##3DMC(a,b) make DMC file later need atleast 4x coverage and methylaton called
	
    ameth = sum(pop_a)
    bmeth = sum(pop_b)

    aunmeth = (len(pop_a) - ameth)
    bunmeth = (len(pop_b) - bmeth)

    if pop_a != pop_b:
    	hstat, pval = s.fisher_exact([[ameth,bmeth],[aunmeth,bunmeth]]) 
	##if sum(pop_a)/float(len(pop_a)) >= .7 and sum(pop_b)/float(len(pop_b)) <= .1:
	if pval <= 0.001:
		print "{0}\t{1}\t{2}\t{3}\t{4}\t{7}\t{5}\t{6}".format(seqid,start,stop,ameth,bmeth,hstat,pval,len(pop_a))
    	##if sum(pop_b)/float(len(pop_b)) >= .7 and sum(pop_a)/float(len(pop_a)) <= .1:
	####	print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(seqid,start,stop,ameth,bmeth,hstat,pval)
    

def main(features,meth_file):
    meth_data = loadintointersect_meth(meth_file)
    freq(features,100,50,meth_data)


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("--genomic", dest="genomic_bed", help="bedfile containing genomic info ie: chr start stops")
    parser.add_option("--meth", dest="meth_file", help="bedfile containing data information ie histone mark mapping")
    (options, _) = parser.parse_args()
    main(options.genomic_bed,options.meth_file)
