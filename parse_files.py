import numpy as np
from collections import defaultdict

def parse_files(input_file, method, all_names_scientific, all_names_other, nodes, merged, ranks, verbose):
	
	def retrieveValidTaxID(taxid, name, rank, method):
		# Normalize taxids and names into one single version of the taxonomic database (one provided as a parameter)
		# 1) First check for valid ranks (provided as variable ranks)
		# 2) Look up for the taxid (nodes.dmp), meaning it's still valid
		# 3) Look up for changes in the taxid (merged.dmp) and return updated version (check if rank is still valid)
		# 4) Look up for scientific name (names.dmp), return only if find unique match
		# 5) Look up for any other taxonomic name (names.dmp), return only if find unique match
 		
		# Verify if rank is valid -> only on profiling (binning is kept to account for the abundance estimation)
		if rank not in ranks.ranks and method=='p': 
			if verbose: print(("Ignored entry [%s] rank [%s] - Rank not expected" % (taxid,rank)))
			return 0
		
		if taxid:
			if taxid in nodes: # Valid taxid
				return taxid
			elif taxid in merged: # Search for updated taxid on merged.dmp
				if verbose: print(("(WARNING) merged taxid [%d] rank [%s] -> [%d] rank [%s]") % (taxid, rank, merged[taxid], nodes[merged[taxid]]['rank']))
				# If found on merged.dmp, verify if rank is valid
				if nodes[merged[taxid]]['rank'] not in ranks.ranks and method=='p':  
					if verbose: print(("Ignored entry [%s] rank [%s] - Rank not expected" % (taxid,rank)))
					return 0
				else:
					return merged[taxid]
		
		if name:
			# If there's no valid taxid, look for name
			if (name,rank) in all_names_scientific:
				# First look for scientific names
				if len(all_names_scientific[(name,rank)])>1: #Name found but with more than one taxid (ambiguous)
					if verbose: print(("Ignored entry [%s] rank [%s] - ambiguous scientific name" % (name,rank)))
					return 0
				else:
					return all_names_scientific[(name,rank)][0] # Unique name found, return taxid
			elif (name,rank) in all_names_other:
				# Check for non-scientific names
				if len(all_names_other[(name,rank)])>1: #Name found but with more than one taxid (ambiguous)
					if verbose: print(("Ignored entry [%s] rank [%s] - ambiguous name" % (name,rank)))
					return 0
				else:
					return all_names_other[(name,rank)][0] # Unique name found, return taxid
			else:	
				if verbose: print(("Ignored entry [%s] rank [%s] - Taxid and name not found" % (name,rank)))
				return 0
		

	def isbioboxes(input_file):
		with open(input_file, 'r') as f: first_line = f.readline()
		return True if first_line.startswith("@") or first_line.startswith("#") else False
	
	def parse_profiling_bioboxes(input_file):
		result = []
		count = defaultdict(lambda: {'total':0,'ignored':0})
		with open(input_file,'r') as f:
			for line in f:
				if line[0]=="@" or line[0]=="#" or line[0]=="\n": continue
				fields = line.rstrip().split('\t')
				taxid = fields[0]
				rank = fields[1]
				name = fields[3].split("|")[-1]
				ab = fields[4]
				
				count[rank]['total']+=1
				taxid = int(taxid) if taxid.isdigit() else 0
				
				valid_taxid = retrieveValidTaxID(taxid, name, rank, "p")
				if valid_taxid: # Add taxid to results
					result.append([1,ranks.getRankID(rank), valid_taxid, float(ab)])
				else:
					count[rank]['ignored']+=1
		return result, count
		
	
	def parse_profiling_tsv(input_file):
		result = []
		count = defaultdict(lambda: {'total':0,'ignored':0})
		with open(input_file,'r') as f:
			for line in f:
				rank, name_taxid, ab = line.rstrip().split('\t')
				
				count[rank]['total']+=1
				taxid = int(name_taxid) if name_taxid.isdigit() else 0
				name = name_taxid if not name_taxid.isdigit() else ""
				
				valid_taxid = retrieveValidTaxID(taxid, name, rank, "p")
				if valid_taxid: # Add taxid to results
					result.append([1,ranks.getRankID(rank), valid_taxid, float(ab)])
				else:
					count[rank]['ignored']+=1
		return result, count
		
	def parse_binning_bioboxes(input_file):
		result = []
		count = {'total':0,'ignored':0}
		with open(input_file,'r') as f:
			for line in f:
				if line[0]=="@" or line[0]=="#" or line[0]=="\n": continue
				fields = line.rstrip().split('\t')
				taxid = int(fields[1])
				len = fields[2]
				count['total']+=1
				valid_taxid = retrieveValidTaxID(taxid, "", "", "b")
				if valid_taxid: # Add taxid to results
					result.append([int(valid_taxid),int(len)])
				else:
					count['ignored']+=1
					if verbose: print(("Ignored entry [%s] - Taxid not found" % (taxid)))

		return result, count

	def parse_binning_tsv(input_file):
		result = []
		count = {'total':0,'ignored':0}
		with open(input_file,'r') as f:
			for line in f:
				fields = line.rstrip().split('\t')
				taxid = int(fields[1])
				len = fields[2]
				count['total']+=1
				valid_taxid = retrieveValidTaxID(taxid, "", "", "b")
				if valid_taxid: # Add taxid to results
					result.append([int(valid_taxid),int(len)])
				else:
					count['ignored']+=1
					if verbose: print(("Ignored entry [%s] - Taxid not found" % (taxid)))
					
		return result, count
	
	def b2p_count(result,nodes):
		# Binning to Profiling
		count = defaultdict(int)
		
		# Sum read lengths for each entry (by taxid)
		count = defaultdict(int)
		for res in result: count[res[0]]+=res[1]
			
		# From the counts (only leafs) reconstruct the taxonomic tree, summing up lengths to higher taxonomic levels
		for taxid,len in list(count.items()):		
			#Start from the parent node because the entry itself (leaf) already contains the count
			txid = nodes[taxid]['parent']
			nlen = len
			while txid!=1: # Until root node
				nlen = count[txid] + nlen # New length account for node accumulated length plus parent length
				count[txid] = nlen # Add new lenght to the parent
				txid = nodes[txid]['parent'] # Set next node
		
		# Generate the profile
		r = []
		for taxid,len in list(count.items()):
			rank = nodes[taxid]['rank']
			# Only keep choosen ranks - They are ignored here but their counts were already summed up to the valid ranks
			if rank in ranks.ranks:
				count[rank]+=1
				r.append([1,ranks.getRankID(rank),taxid,float(len)])
		
		return r, count
	
	if method=='db' or method=='p':
		if isbioboxes(input_file):
			print(" - %s (BioBoxes)" % input_file)
			r, count = parse_profiling_bioboxes(input_file)
		else:
			print(" - %s (tsv)" % input_file)
			r, count = parse_profiling_tsv(input_file)
		tot = 0
		ign = 0
		for rank in ranks.ranks: 
			tot+=count[rank]['total']
			ign+=count[rank]['ignored']
			print(("\t%s - %d entries (%d ignored)") % (rank, count[rank]['total']-count[rank]['ignored'],count[rank]['ignored']))
		print(("\tTotal - %d taxons (%d ignored)") % (tot,ign))
	elif method=='b':
		if isbioboxes(input_file):
			print(" - %s (BioBoxes)" % input_file)
			br, count = parse_binning_bioboxes(input_file)
			r, count2 = b2p_count(br,nodes)
		else:
			print(" - %s (tsv)" % input_file)
			br, count = parse_binning_tsv(input_file)
			r, count2  = b2p_count(br,nodes)
		print(("\t%d lines (%d ignored)") % (count['total'],count['ignored']))
		tot = 0
		for rank in ranks.ranks: 
			tot+=count2[rank]
			print(("\t%s - %d entries") % (rank, count2[rank]))
		print(("\tTotal - %d taxons") % (tot))

	#r = [Presence,RankID,TaxID,Val]
	return np.array(r)
