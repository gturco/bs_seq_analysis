from intersection import Intersecter, Feature
from flatfeature import Bed
from collections import defaultdict
from random import shuffle
import re
from pyfasta import Fasta

def loadintointersect_meth(meth_file):
    query_list_pos = {}
    query_list_neg = {}
    for line in open(meth_file):
        seqid, start, strand = line.strip().split("\t")[0:3]
        if seqid[0] != "c": continue
        seqid = seqid.strip("chromosome_") 
        seqid = seqid.strip()
	if strand == "+":
    ### ADD one because bed adds one too number
            if seqid not in list(query_list_pos.keys()):
         		query_list_pos[seqid] = Intersecter()
            query_list_pos[seqid].add_interval(Feature(int(start), int(start),name=strand))
        elif strand == "-":
            if seqid not in list(query_list_neg.keys()):
         		query_list_neg[seqid] = Intersecter()
            query_list_neg[seqid].add_interval(Feature(int(start), int(start),name=strand)) 
    return query_list_pos, query_list_neg



def loadintointersect(bed_file):
    query_list_pos = {}
    query_list_neg = {}
    feature_list = Bed(bed_file)
    for feature in feature_list:
    ##    if float(feature['accn']) < .4: continue
        if feature["strand"] == "+":
    ### ADD one because bed adds one too number
            if feature['seqid'] not in list(query_list_pos):
         		query_list_pos[feature['seqid']] = Intersecter()
            query_list_pos[feature['seqid']].add_interval(Feature(int(feature['start']-1), int(feature['start']-1),name=feature['strand']))
        elif feature["strand"] == "-":
            if feature['seqid'] not in list(query_list_neg):
         		query_list_neg[feature['seqid']] = Intersecter()
            query_list_neg[feature['seqid']].add_interval(Feature(int(feature['start']-1), int(feature['start']-1),name=feature['strand'])) 
    return query_list_pos, query_list_neg

def get_matchs(query_list,seqid,searchrange,fasta,genomic_location,rc):
    """ searches through a arange counting the number of meth sites in each window and returning the location
    and number of sites there"""
    positions = []
    sites = []
    window_size = 20
    pos = genomic_location
    for region_start in searchrange[::window_size]:
        region_end = region_start + window_size
    #cs = len(cs_list[seqid].find(region_start,region_end))
    #cs = get_cs(fasta,region_start,region_end,rc) 
        matches = query_list[seqid].find(region_start, region_end)
        sites.append(len(matches))
        #if len(matches) == 0 :
	#	sites.append(0)
	#elif cs == 0:
 	#	print "ERROR",len(matches),region_start,region_end,seqid,cs
	
	
	positions.append(pos)
	pos += window_size

    if rc:
	return zip(positions,sites[::-1])
    else: 
	return zip(positions,sites)

def get_exonmeth(feature,query_list):
    exon_meth = []
    for e in feature['locs']:
        for i in range(e[0],e[1]+1):
	    matches = query_list[feature['seqid']].find(i, i)
        #c = len(cs_list[feature['seqid']].find(i, i))
            exon_meth.append(len(matches))
    	#cs.append(c)
        ##seq += fasta[e[0]:e[1]+1]
    ###print len(exon_meth), sum(exon_meth)
    return exon_meth


def get_genebody(query_list,feature,fasta,rc,rand):
    sites = []
    pos = []
    exon_meth = get_exonmeth(feature,query_list)
    ebin = len(exon_meth)/100.00
    index = 0
    ebin_n = 0
    #### make a randmom about add some
    dec = (ebin - int(ebin))
    add_one = int(dec * 100)
    if len(rand) > 0:
	randomize = rand[feature["accn"]]
    else:
    	randomize = add_one * [False] + (100-add_one) * [True]
    	shuffle(randomize)


    for n in range(0,100):
        if randomize[n]:
            ebin_n += int(ebin)
            m = sum(exon_meth[index:int(ebin_n)])
        else:
            ebin_n += int(ebin) + 1
            m = sum(exon_meth[index:int(ebin_n)])
 	##region = seq[index:int(ebin_n)]
    #if rc : cs = re.findall("G",region)
	#else: cs = re.findall("C",region)
	#if m == 0 :
	#	sites.append(0)
	#elif cs == 0:
	#	print "ERROR",m,index,ebin_n,feature["accn"],ebin,exon_meth[index:ebin_n],seq
	#else:        

	sites.append(float(m))
	pos.append(10*n)
        index = ebin_n
    ##print sum(sites), sum(exon_meth)   
    if rc :
    	return zip(pos,sites[::-1]), randomize
    else: 
	return zip(pos,sites), randomize

def main(feature_bed,query_list_pos,query_list_neg,fasta_file,mtype,rand):
    features = Bed(feature_bed)
    fasta = Fasta(fasta_file)
    All_sites = defaultdict(list)
    r = {}
    cgene = {}
    for feature in features:
	rc = feature["strand"] == "-"
	if feature["strand"] == "+":
        	TSS_region = range(int(feature['locs'][0][0]) - 2000,int(feature['locs'][0][0]))
        	TTS_region = range(int(feature['locs'][-1][1]),int(feature['locs'][-1][1]) + 2000)
        	TSS_sites = get_matchs(query_list_pos,feature['seqid'],TSS_region,fasta["chromosome_" + feature["seqid"]],-2000,rc)
        	TE_sites = get_matchs(query_list_pos,feature['seqid'],TTS_region,fasta["chromosome_" + feature['seqid']],1100,rc)
        	gene_body, rebin = get_genebody(query_list_pos,feature,fasta["chromosome_" + feature["seqid"]],rc,rand)
		r[feature["accn"]] = rebin	
		cgene[feature["accn"]] =gene_body	
                
         #       [All_sites[str(region)].append(freq) for region,freq in TSS_sites]
        # 	[All_sites[str(region)].append(freq) for region,freq in TE_sites]
		[All_sites[feature["accn"]].append((region,freq)) for region,freq in TSS_sites]
		[All_sites[feature["accn"]].append((region,freq)) for region,freq in TE_sites]




	if feature["strand"] == "-":
		TTS_region = range(int(feature['locs'][0][0]) - 2000,int(feature['locs'][0][0]))
        	TSS_region = range(int(feature['locs'][-1][1]),int(feature['locs'][-1][1]) + 2000)
        	TSS_sites = get_matchs(query_list_neg,feature['seqid'],TSS_region,fasta["chromosome_" + feature["seqid"]],-2000,rc)
        	TE_sites = get_matchs(query_list_neg,feature['seqid'],TTS_region,fasta["chromosome_" + feature['seqid']],1100,rc)
        	
		###RV complent
		gene_body,rebin = get_genebody(query_list_neg,feature,fasta["chromosome_" + feature["seqid"]],rc,rand)
		r[feature["accn"]] = rebin	
		cgene[feature["accn"]] = gene_body	
        	

        	##[All_sites[str(region)].append(freq) for region,freq in TSS_sites]
        	##[All_sites[str(region)].append(freq) for region,freq in TE_sites]
         	[All_sites[feature["accn"]].append((region,freq)) for region,freq in TSS_sites]
		[All_sites[feature["accn"]].append((region,freq)) for region,freq in TE_sites]


    return All_sites,r,cgene


def get_cs(fasta,start,end,rc):
	seq = fasta[start:end+1]
	Cs = re.findall("C",seq)
	if rc:
		Cs = re.findall("G",seq)
        return len(Cs)

def maketable(dic,cs,cgene,mgene,mtype):
    #outfile = open(filename, "wb")
##    header = list(dic.keys())
##    values = []
##    #outfile.write("pos\tfreq\n")
##    for key in dic.keys():
##        try:
##		value = sum(dic[key])/float(sum(cs[key]))
##        except ZeroDivisionError:
##		if sum(dic[key]) == 0:
##			value = 0
##		else: value = "NA"
##	values.append(value)
##        line = "{0}\t{1}\t{2}".format(key,value,mtype)
###        line = "{0}\t{1}\t{2}".format(key,sum(cs[key]),mtype)
##        
##        #outfile.write(line)
##    	print line
    d = defaultdict(list)
    for gene_name in cgene.keys():
        for i,bin in enumerate(cgene[gene_name]):
        	region, freq = bin
        	if freq != 0 :
        		d[region].append(mgene[gene_name][i][1]/float(freq))
        for i,bin in enumerate(cs[gene_name]):
        	cregion, cfreq = bin
        	if cfreq != 0 :
        		d[cregion].append(dic[gene_name][i][1]/float(cfreq))

    for cregion in d.keys():
 ###       if len(d[cregion]) < 3: continue
	value = sum(d[cregion])/float(len(d[cregion]))
        line = "{0}\t{1}\t{2}".format(cregion,value,mtype)
   	print line
 
 #outfile.close()

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("-m", dest="meth_file", help=" containing ALL Cs")
    parser.add_option("-f", dest="methfile", help="bedfile containing sig hits islandfiltered.bed")
    parser.add_option("--ref", dest="ref", help="ref genome bed (sorg_v2)")
    parser.add_option("--fasta", dest="fasta", help="ref genome fasta (sorg_v2)") 
    parser.add_option("--mtype", dest="mtype", help="methylation type")
    (options, _) = parser.parse_args()


    query_list_pos,query_list_neg = loadintointersect_meth(options.meth_file)
    y = []
    cs,rand,cgene = main(options.ref, query_list_pos, query_list_neg , options.fasta  ,options.mtype, y)
    del query_list_neg, query_list_pos
    query_list_pos,query_list_neg = loadintointersect_meth(options.methfile)
    meth,rand,mgene = main(options.ref, query_list_pos, query_list_neg , options.fasta  ,options.mtype, rand)
    maketable(meth,cs,cgene,mgene,options.mtype)
