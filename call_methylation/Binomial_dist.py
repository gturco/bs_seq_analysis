from scipy.stats import binom_test


def get_fdr(p_val_distribution,uncorrected,alpha):

    q = sum(1 for x in p_val_distribution if x < uncorrected) \
            * 1./len(p_val_distribution)
    return q

def filter_methylhits(meth_file,p):
    lines = []
    distribution = []
    print "seqid\tstart\tstrand\tmeth_type\tfraction\ttotal\tpval"
    for line in open(meth_file):
        seqid,start,strand,meth_type,fraction,total = line.strip().split("\t")
        total = int(total)
        methcount = float(eval(fraction)) * int(total)
        if total < 4 or methcount == 0: continue
        pval = binom_test(methcount,int(total),p)
        distribution.append(pval)
        print line.strip() + "\t{0}".format(pval)
        lines.append(line)

def main(meth_file,p):
    filter_methylhits(meth_file,p)
    ####run_r_script

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("--meth", dest="meth", help="methfile from rolf")
    parser.add_option("--pval", dest="pval",type='float', help="spike in conrol error rate")
    (options, _) = parser.parse_args()
    
    main(options.meth,options.pval)
