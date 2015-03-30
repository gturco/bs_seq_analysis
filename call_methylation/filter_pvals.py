
pvalfile = "all_nonvascular.mr.meth.CHH.hypo.1.pvalcorrected"


for line in open(pvalfile):
	if line[:2] == "se":
		print line.strip()
	else: 	
		reads = line.split("\t")[5]
		if int(reads) < 4 : continue
		print line.strip()

