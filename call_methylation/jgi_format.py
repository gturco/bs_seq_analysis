def jgi_format(methyl_file,out_file):
	out = open(out_file, "wb")
	for line in open(methyl_file):
		### if bismark header skip
        	if line[0] == "s": continue
		seqid, start, strand, mytpe, score = line.strip("\n").split("\t")[:5]
		new_line = "{0}\t{1}\t{2}\t{3}\t1\t{4}\n".format(seqid.strip("chromosome_"),start,int(start)+1, score,strand)
		out.write(new_line)
	out.close()

###3jgi_format("all_nonvascular/all_nonvascular.mr.meth.CHH.hypo.pvalcorrected","all.nonvascular.CHH.bed") 
###jgi_format("all_root_CHH/all.root.mr.meth.CHH.hypo.pvalcorrected","all.root.CHH.bed") 
###jgi_format("all_shoot_CHG/all.shoot.mr.meth.CHG.hypo.pvalcorrected","all.shoot.CHG.bed") 
###jgi_format("all_shootCHH/all.shoot.mr.meth.CHH.hypo.pvalcorrected","all.shoot.CHH.bed") 
###jgi_format("nonvas_CG/all_nonvascular.mr.meth.CG.hypo.pvalcorrected","all.nonvas.CG.bed") 
jgi_format("vascular_CHG/vascular.mr.meth.CHG.hypo.pvalcorrected","all.vascular.CHG.bed") 
###jgi_format("vascular_CHH/vascular.mr.meth.CHH.hypo.pvalcorrected","all.vascular.CHH.bed") 
jgi_format("all_vas_CG/vascular.mr.meth.CG.hypo","all.vascular.CG.bed") 


