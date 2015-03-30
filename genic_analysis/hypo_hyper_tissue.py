from collections import defaultdict

def main(meth_file,c_type):
	d = defaultdict(list)
	vcs = []
	nvc = []	
	for i,line in enumerate(open(meth_file)):
		if i % float(20) == 0:
			change = float(len(nvc)) - len(vcs)
			per = change/float(20)
			###print i ,change , vcs, nvc
			##if change > 0 or change < 0:
			###	print i,change, nvc, vcs
			d[per].append(i)
			vcs = []
			nvc = []
			print per
		else:
			
			args = line.strip().split("\t")
			if float(args[4]) > 0:
				nvc.append(1)
			if float(args[10]) > 0:
				vcs.append(1)
		
	####for key in d.keys():
	####	per = len(d[key])/float(i)
	####	print "{0}\t{1}\t{2}".format(key,per,c_type)


##main("../non_vas_CG_master.txt","CG")
##main("../non_vas_CHG_master.txt","CHG")
###main("../all.non_vascular.CHH.master.txt","CHH")


###main("../DMR/root_shoot_CG_master.txt","CG")
###main("../DMR/root_shoot_CHG_master.txt","CHG")
###main("../DMR/all.root_shoot.CHH.master.txt","CHH")

main("../DMR/root_vas_CG_master.txt","CG")
main("../DMR/root_vas_CHG_master.txt","CHG")
main("../DMR/all.root_vascular.CHH.master.txt","CHH")

###main("../DMR/root_non_CG_master.txt","CG")
###main("../DMR/root_non_CHG_master.txt","CHG")
###main("../DMR/all.root_nonvascular.CHH.master.txt","CHH")

