from collections import defaultdict

def parse_tax(names_file,nodes_file,merged_file,ranks):
	all_names_scientific = defaultdict(list)
	all_names_other = defaultdict(list)
	nodes = {}
	merged = {}
	
	for l in open(nodes_file,'r'):
		taxid, parent_taxid, rank, _ = l.split('\t|\t',3)
		nodes[int(taxid)] = {'parent':int(parent_taxid),'rank':rank,'name':''}
	
	for l in open(merged_file,'r'):
		old_taxid, new_taxid, _ = l.rstrip().split('\t|',2)
		merged[int(old_taxid)] = int(new_taxid)
	
	for l in open(names_file,'r'):
		fields = l.split('\t|\t')
		taxid = int(fields[0])
		rank = nodes[taxid]['rank']
		
		# Only parse selected ranks
		if rank in ranks.ranks:
			name = fields[1]
			nc = fields[3].replace('\t|\n','')
		
			if nc=="scientific name":
				if (name,rank) in all_names_scientific: print("parse_tax WARNING repeated SCIENTIFIC (name,rank)", rank, name, taxid, all_names_scientific[(name,rank)])
				all_names_scientific[(name,rank)].append(taxid)
			else:
				if (name,rank) in all_names_other: print("parse_tax WARNING repeated OTHER (name,rank)", rank, name, taxid, all_names_other[(name,rank)])
				all_names_other[(name,rank)].append(taxid)

			# Set scientific name to taxid (or any other in case is still empty)
			if nc=="scientific name" or not nodes[taxid]['name']:
				nodes[taxid]['name'] = fields[1]

	return all_names_scientific, all_names_other, nodes, merged
