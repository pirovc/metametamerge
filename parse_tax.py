from collections import defaultdict

def parse_tax(names_file,nodes_file,ranks):
	all_names_scientific = defaultdict(list)
	all_names_other = defaultdict(list)
	nodes = {}
	
	for l in open(nodes_file,'r'):
		fields = l.split('\t|\t')
		nodes[int(fields[0])] = {'parent':int(fields[1]),'rank':fields[2],'name':''}
	
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

	return all_names_scientific, all_names_other, nodes
