python Binomial_dist.py --meth "all_nonvascular.mr.meth.CHH.hypo.1" --pval 0.001 > "all_nonvascular.mr.meth.CHH.hypo.1.pval"
R CMD BATCH meth_fdr.R
python filter_pvals.py
