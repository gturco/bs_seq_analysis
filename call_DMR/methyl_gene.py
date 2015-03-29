import csv
from intersection import Intersecter, Feature
from collections import defaultdict
import sys
from flatfeature import Bed

def index_data(csv_path):
    data = csv.DictReader(open(csv_path, 'rb'), delimiter='\t', fieldnames = ['seqid','start','stop','name',"number","strand"])
    return data


def group_by_seqids(dfile):
    """group lines by seqid {seqid:[{line1},{line2}]}"""
    #### NOT ALWAYS INT
    gseqids = defaultdict(list)
    for line in dfile:
        try:
            seqid = int(line['seqid'].strip("chromosome_"))
        except ValueError: seqid = 11
        gseqids[seqid].append(line)
    return gseqids

def insert_queries(quey_file):
    """creates an interval list containg all the peaks"""
    d = index_data(quey_file)
    gseqids = group_by_seqids(d)
    seqid_numb = len(gseqids.keys())
    query_list = []
    for n in range(0,seqid_numb):
        query_list.append(Intersecter())
    for seqid in range(0,seqid_numb):
        #inserts query list into intervals
        #print gseqids[4]
        [query_list[seqid].add_interval(Feature(int(query_info['start']),
        int(query_info['stop']),name=query_info['number'],strand=int(query_info['strand']+"1"))) for query_info in gseqids[seqid+1]]
    return query_list


def get_genespace(bed,locs,gstart,gend):
    feats = bed.get_features_in_region(str(locs['seqid']), locs['start'], locs['end'])
    feat_starts = [int(f['start']) for f in feats if not (f['start'] == locs['start'] and f['end'] == locs['end'])]
    feat_ends = [int(f["end"]) for f in feats if not (f['start'] == locs['start'] and f['end'] == locs['end'])]
    feat_starts.sort()
    feat_ends.sort()
    if len(feats) > 1:
        if feat_ends[0] < int(locs['start']):
             gene_spacestart = feat_ends[0]
        else:
            gene_spacestart = gstart
        if feat_starts[-1] > int(locs['end']):
            gene_spacend = feat_starts[-1]
        else:
            gene_spacend = gend
    else:
        gene_spacend = gend
        gene_spacestart = gstart

    return gene_spacestart, gene_spacend

def amt_overlap(matches,strand,pad_start,pad_stop):
    """gives amount of overlap"""
    overlap = []
    real_m = []
    p_m = []
    for match in matches:
        if strand != match.strand:
            overlap.append(0)
            continue
        ed = match.stop >= pad_stop
        st = pad_start >= match.start
        if st and ed:
             overlap.append(pad_stop - pad_start)
             real_m.append(match.name)
        elif not st and not ed:
            overlap.append(match.stop - match.start)
        elif st:
            overlap.append(match.stop - pad_start)
        elif ed:
            overlap.append(pad_stop - match.start)
    if len(overlap) == 0:
        return 0
    else:
        return sum(overlap)/float(pad_stop-pad_start)

def find_intersections(padding,interval_list,genelist,outfh):
    "takes a list of gene postions and compares it to chipseq data \
            and returns value of peaks that overlap with CNS" 
    """ADD columns number_same, number_diff, match_same, match_diff"""
    out = open(outfh,'wb')
    header ='seqid,start,stop,name,peak_score\n'
    for gene in genelist:
        ####get_genespace
        gene = genelist.row_to_dict(gene)
        pad_start = max(0,int(gene['start']) - padding)
        pad_stop = int(gene['end']) + padding
        pad_start,pad_stop = get_genespace(genelist,gene,pad_start,pad_stop)
        three_prom = interval_list[int(gene['seqid'])-1].find((pad_start),int(gene['start']))
        five_prom = interval_list[int(gene['seqid'])-1].find(( int(gene['end'])),(pad_stop))
        gene_body = interval_list[int(gene['seqid'])-1].find((int(gene['start'])),(int(gene['end'])))
        
        strand = int(gene['strand'] + "1")
        three_prom_p = amt_overlap(three_prom,strand,pad_start,int(gene['start']))
        five_prom_p = amt_overlap( five_prom, strand, int(gene['end']),pad_stop)
        gene_body_p = amt_overlap(gene_body,strand,int(gene['start']),int(gene['end']))

        gene_name = gene["accn"]
       # if gene_name == "Sb04g003710":
       #     print gene
       #     print interval_list[0].find(0,100000000)
       #     print interval_list[0].find(3577840,3577841)
       #     print three_prom
       #     print gene_body
       #     print five_prom
        #three_prom = [i for i in three_prom if i.name == gene['strand']]
        #five_prom = [i for i in five_prom if i.name == gene['strand']]
        #gene_body = [i for i in gene_body if i.name == gene['strand']]
        
       

        if len(three_prom) > 0:
            l = "{0}\t3_prom\t{1}\t{2}\t{3}\n".format(gene_name,three_prom_p,three_prom_p *len(three_prom),sum(int(sig.name) for sig in three_prom))
            out.write(l)
        if len(five_prom) > 0:
            l = "{0}\t5_prom\t{1}\t{2}\t{3}\n".format(gene_name,five_prom_p,five_prom_p * len(five_prom),sum(int(sig.name) for sig in five_prom))
            out.write(l)
        if len(gene_body) > 0:
            l = "{0}\tgene_body\t{1}\t{2}\t{3}\n".format(gene_name,gene_body_p,gene_body_p * len(gene_body),sum(int(sig.name) for sig in gene_body))
            out.write(l)


genelist = Bed("sorg.bed")
interval_list = insert_queries("DMR_NONVAS_CG_HYPO")

find_intersections(3000,interval_list,genelist,"DMR_nonvas.genes")
