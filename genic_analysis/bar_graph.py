from intersection import Intersecter, Feature
from flatfeature import Bed
from collections import defaultdict



def loadintointersect_meth(meth_file):
    query_list_pos = {}
    query_list_neg = {}
    for line in open(meth_file):
        seqid, start, strand = line.strip().split("\t")[0:3]
        if seqid[0] != "c": continue
        seqid = seqid.strip("chromosome_") 
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




def main_gene(feature_file, query_list_pos,query_list_neg):
    cds = []
    three_all = []
    five_all = []
    feature_bed = Bed(feature_file)
    for feature in feature_bed:
        exon_meth = []
        for e in feature['locs']:
            for i in range(e[0],e[1]+1):
                if feature["strand"] == "+": 
                    matches = query_list_pos[feature['seqid']].find(i, i)
                else:
                    matches = query_list_neg[feature['seqid']].find(i, i)
 
                exon_meth.append(len(matches))
        cds.append(sum(exon_meth))
        if feature["strand"] == "+":
            five_prime = query_list_pos[feature['seqid']].find(int(feature['locs'][0][0])-300, int(feature['locs'][0][0]))
            three_prime = query_list_pos[feature['seqid']].find(int(feature['locs'][-1][1]),int(feature['locs'][-1][1]) + 300)
        elif feature["strand"] == "-":
            three_prime = query_list_pos[feature['seqid']].find(int(feature['locs'][0][0])-300, int(feature['locs'][0][0]))
            five_prime = query_list_pos[feature['seqid']].find(int(feature['locs'][-1][1]),int(feature['locs'][-1][1]) + 300)
        three_all.append(len(three_prime))
        five_all.append(len(five_prime))
    return cds, three_all, five_all



def main(feature_file, query_list_pos, query_list_neg):
    freq = []
    for feature in open(feature_file):
        seqid,start,end,name,strand = feature.strip().split("\t")
        start, end = int(start), int(end)
	if strand == "+":
            m = query_list_pos[seqid].find(start,end)
        else:
            m = query_list_neg[seqid].find(start,end)
        freq.append(len(m))
    return freq

def maketable(feature_freq,cov_freq,ttype,mtype):
    #outfile = open(filename, "wb")
    values = []
    #outfile.write("pos\tfreq\n")
    for i,freq in enumerate(feature_freq):
        try:
	    value = (freq)/float(cov_freq[i])
        except ZeroDivisionError: continue
        values.append(value)
    line = "{0}\t{1}\t{2}".format(sum(values)/float(len(values)),ttype,mtype)
    print line


def run_all(meth_file,cov_file,sorg_bed,mirna_bed,trans_bed,dups_bed,mtype):

    query_list_pos, query_list_neg = loadintointersect_meth(meth_file)
    mirna = main(mirna_bed,query_list_pos,query_list_neg)
    trans = main(trans_bed,query_list_pos,query_list_neg)
    dups = main(dups_bed,query_list_pos,query_list_neg)
    gene_body, three_prime, five_prime = main_gene(sorg_bed,query_list_pos,query_list_neg)





    del query_list_neg, query_list_pos
    query_list_pos, query_list_neg = loadintointersect_meth(cov_file)
    mirna_cov = main(mirna_bed,query_list_pos,query_list_neg)
    trans_cov = main(trans_bed,query_list_pos,query_list_neg)
    dups_cov = main(dups_bed,query_list_pos,query_list_neg)
    gene_body_cov, three_prime_cov, five_prime_cov = main_gene(sorg_bed,query_list_pos,query_list_neg)
    
    
    
    
    maketable(gene_body,gene_body_cov,"CDs",mtype)
    maketable(three_prime,three_prime_cov,"3 promoter",mtype)
    maketable(five_prime,five_prime_cov,"5 promoter",mtype)
    maketable(mirna,mirna_cov,"miRNA",mtype)
    maketable(trans, trans_cov, "TE",mtype)
    maketable(dups,dups_cov,"tandem duplicates",mtype)
    
###run_all("/share/brady-archive1/gturco/sorg_XSEDE/figure2/breakup/all_vas_CG/vascular.mr.meth.CG.hypo.pvalcorrected.filtered","/share/brady-archive1/gturco/sorg_XSEDE/circos/vas_CG_cov_meth","sorg.bed","miRNA.bed","trans.bed","dup.bed","CG")
###run_all("/share/brady-archive1/gturco/sorg_XSEDE/figure2/breakup/vascular_CHG/vascular.mr.meth.CHG.hypo.pvalcorrected.filtered", "/share/brady-archive1/gturco/sorg_XSEDE/figure2/TSS/vas_CHG_cov_meth","sorg.bed","miRNA.bed","trans.bed","dup.bed","CHG")

####run_all("/share/brady-archive1/gturco/sorg_XSEDE/figure2/breakup/vascular_CHH/vascular.mr.meth.CHH.hypo.pvalcorrected.filtered", "/share/brady-archive1/gturco/sorg_XSEDE/figure2/TSS/vas_CHH_cov_meth","sorg.bed","miRNA.bed","trans.bed","dup.bed","CHH")


###################

###run_all("/share/brady-archive1/gturco/sorg_XSEDE/figure2/breakup/nonvas_CG/all_nonvascular.mr.meth.CG.hypo.pvalcorrected.filtered","/share/brady-archive1/gturco/sorg_XSEDE/figure2/TSS/nonvas_CG_cov_meth","sorg.bed","miRNA.bed","trans.bed","dup.bed","CG")
###run_all("/share/brady-archive1/gturco/sorg_XSEDE/figure2/pvals/nonvascular_CHG/all.nonvascular.mr.meth.CHG.hypo.pvalcorrected.filtered", "/share/brady-archive1/gturco/sorg_XSEDE/figure2/TSS/nonvas_CHG_cov_meth","sorg.bed","miRNA.bed","trans.bed","dup.bed","CHG")

###run_all("/share/brady-archive1/gturco/sorg_XSEDE/figure2/breakup/all_nonvascular/all_nonvascular.mr.meth.CHH.hypo.pvalcorrected.filtered", "/share/brady-archive1/gturco/sorg_XSEDE/figure2/TSS/nonvas_CHH_cov_meth","sorg.bed","miRNA.bed","trans.bed","dup.bed","CHH")

###################
### ROOt
###run_all("/share/brady-archive1/gturco/sorg_XSEDE/figure2/pvals/root_CG/all.root.mr.meth.CG.hypo.pvalcorrected.filtered","/share/brady-archive1/gturco/sorg_XSEDE/figure2/TSS/root_CG_cov_meth","sorg.bed","miRNA.bed","trans.bed","dup.bed","CG")
###run_all("/share/brady-archive1/gturco/sorg_XSEDE/figure2/pvals/root_CHG/all.root.mr.meth.CHG.hypo.pvalcorrected.filtered", "/share/brady-archive1/gturco/sorg_XSEDE/figure2/TSS/root_CHG_cov_meth","sorg.bed","miRNA.bed","trans.bed","dup.bed","CHG")

###run_all("/share/brady-archive1/gturco/sorg_XSEDE/figure2/breakup/all_root_CHH/all.root.mr.meth.CHH.hypo.pvalcorrected.filtered", "/share/brady-archive1/gturco/sorg_XSEDE/figure2/TSS/root_CHH_cov_meth","sorg.bed","miRNA.bed","trans.bed","dup.bed","CHH")


###################

run_all("/share/brady-archive1/gturco/sorg_XSEDE/figure2//pvals/shoot_CG/all.shoot.mr.meth.CG.hypo.pvalcorrected.filtered","/share/brady-archive1/gturco/sorg_XSEDE/figure2/TSS/shoot_CG_cov_meth","sorg.bed","miRNA.bed","trans.bed","dup.bed","CG")
run_all("/share/brady-archive1/gturco/sorg_XSEDE/figure2/breakup/all_shoot_CHG/all.shoot.mr.meth.CHG.hypo.pvalcorrected.filtered", "/share/brady-archive1/gturco/sorg_XSEDE/figure2/TSS/shoot_CHG_cov_meth","sorg.bed","miRNA.bed","trans.bed","dup.bed","CHG")

###run_all("/share/brady-archive1/gturco/sorg_XSEDE/figure2/breakup/all_nonvascular/all_nonvascular.mr.meth.CHH.hypo.pvalcorrected.filtered", "/share/brady-archive1/gturco/sorg_XSEDE/figure2/TSS/nonvas_CHH_cov_meth","sorg.bed","miRNA.bed","trans.bed","dup.bed","CHH")





