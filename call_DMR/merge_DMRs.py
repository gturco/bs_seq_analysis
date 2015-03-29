def main(dmr_file):
    dmrs = read_DMR(dmr_file)
    mdmr = merged_dmrs(dmrs)
    for dmr in mdmr:
        print dmr


def read_DMR(dmr_file):
    dmr_data = []
    for line in open(dmr_file):
        seqid,start,end,a,b,total = map(int,line.strip().split("\t")[:5])
        dmr_data.append(seqid,start,end,a,b,total)
    return dmr_data


def merge_dmrs(dmr_data):

    while overlap == True:
        dmrs = dmr_data
        merged_dmrs = []
        overlap = False
        oseqid = 1
        oend = 1
        for dmr in dmrs:
            seqid,start,end,a,b,total = dmr
            if oseqid == seqid and start <= oend:
                merged_dmr = 
                merged_dmr.append(merged_dmr)
                overlap = True
            else:
                dmrs.append(dmr)
            oseqid = seqid
            oend = end
            dmr_data = merged_dmr

    return merged_dmrs

main(dmr_file)
