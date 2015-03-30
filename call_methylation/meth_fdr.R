###R
###meth_fdr.R
##source("http://bioconductor.org/biocLite.R")
###biocLite("multtest")
library(multtest)


fdr_mcount = function(mcount_file,outfile)
{
    mcount = read.table(mcount_file,header=TRUE, sep="\t")

    fdr_v <-  mt.rawp2adjp(mcount$pval,"BH")

    mcount_adj <- mcount[fdr_v$index,]
    mcount_fdr <- cbind(mcount_adj,fdr_v$adjp)

    mcount_0.01_fdr = mcount_fdr[mcount_fdr$BH <= 0.01,]
    write.table(mcount_0.01_fdr[,0:6], outfile ,quote=FALSE, sep="\t",row.names=FALSE)
    
    }


fdr_mcount("all_nonvascular.mr.meth.CHH.hypo.1.pval","all_nonvascular.mr.meth.CHH.hypo.1.pvalcorrected")







