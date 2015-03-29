def jgi_format(methyl_file,out_file):
	out = open(out_file, "wb")
	for line in open(methyl_file):
		### if bismark header skip
        	if line[0] == "s": continue
		seqid, start, strand, mytpe, score,total = line.strip("\n").split("\t")
		if int(total) > 3 and int(total) < 21: 
		####if float(score) > .40:
			out.write(line)
	out.close()

jgi_format("all_nonvascular.mr.meth.CHH.hypo","all_nonvascular.mr.meth.CHH.hypo.resample") 


