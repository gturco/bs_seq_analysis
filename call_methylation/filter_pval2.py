
pvalfile = "all.root.mr.meth.CHG.hypo.pval"


for line in open(pvalfile):
	if line[:2] == "se":
		print line.strip()
	else: 	
		reads = line.split("\t")[5]
		if int(reads) < 4 : continue
		print line.strip()

