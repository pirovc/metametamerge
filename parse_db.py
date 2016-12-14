from collections import defaultdict

def parse_db(input_files,all_names_scientific, all_names_other):
	dbs = defaultdict(int)
	for file in input_files:	
		for line in open(file,'r'):
			rank, name, seq_len = line.rstrip().split('\t')
			if (name,rank) in all_names_scientific:
				if len(all_names_scientific[(name,rank)])>1: print("parse_db WARNING ambiguous name in all_names_scientificy", rank, name)
				#all_names_scientific[(name,rank)][0] - first taxid (can be more than one if repeated entry)
				dbs[all_names_scientific[(name,rank)][0]]+=1 # taxid index
			elif (name,rank) in all_names_other: 
				if len(all_names_other[(name,rank)])>1: print("parse_db WARNING ambiguous name in all_names_other", rank, name)		
				#all_names_other[(name,rank)][0] - first taxid (can be more than one if repeated entry)
				dbs[all_names_other[(name,rank)][0]]+=1 # taxid index
			else: 
				print("parse_db ERROR ignoring name not found in all_names_scientific and all_names_other", rank, name)
	
	return dbs
