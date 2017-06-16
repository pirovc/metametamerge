import numpy as np
from collections import defaultdict

def parse_files(input_file, method, all_names_scientific, all_names_other, nodes, merged, ranks, verbose):
	
	def retrieveValidTaxID(taxid, name=None, rank=None):
		# Normalize taxids into one single version of the taxonomic database (one provided as a parameter)
		# 1) Look up for the taxid (nodes.dmp), meaning it's still valid
		# 2) Look up for changes in the taxid (merged.dmp) and return updated version (check if rank is still valid)
		# 3) Look up for scientific name (names.dmp), return only if find unique match
		# 4) Look up for any other taxonomic name (names.dmp), return only if find unique match
		if taxid:
			if taxid in nodes: # Valid taxid
				return taxid
			elif taxid in merged: # Search for updated taxid on merged.dmp
				if verbose: print(("(WARNING) merged taxid [%d] -> [%d] rank [%s]") % (taxid, merged[taxid], nodes[merged[taxid]]['rank']))
				return merged[taxid]
		
		# If there's no valid taxid, look for name
		if name and rank:
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

		if verbose: print(("Ignored entry [%s] rank [%s] - Taxid and name not found" % (taxid if taxid else name,rank)))
		return 0

	def parse_profiling(input_file, isbiob):
		result = []
		count = defaultdict(lambda: {'total':0,'ignored':0})
		with open(input_file,'r') as f:
			for line in f:
				if isbiob:
					if line[0]=="@" or line[0]=="#" or line[0]=="\n": continue
					fields = line.rstrip().split('\t')
					taxid = int(fields[0]) if fields[0].isdigit() else 0 
					rank = fields[1]
					name = fields[3].split("|")[-1]
					ab = fields[4]
				else:
					rank, name_taxid, ab = line.rstrip().split('\t')
					taxid = int(name_taxid) if name_taxid.isdigit() else 0
					name = name_taxid if not name_taxid.isdigit() else ""
				
				# Verify if rank is valid -> only on profiling (binning is kept to account for the abundance estimation)
				if rank not in ranks.ranks: 
					if verbose: print(("Ignored entry [%s] rank [%s] - Rank not expected" % (taxid if taxid else name,rank)))
					count[rank]['ignored']+=1
				else:
					valid_taxid = retrieveValidTaxID(taxid, name, rank)
					#If taxid changed, re-check rank
					if valid_taxid:
						if valid_taxid!=taxid and nodes[valid_taxid]['rank'] not in ranks.ranks:
							if verbose: print(("Ignored entry [%s] rank [%s] - Rank not expected" % (taxid if taxid else name,rank)))
							count[rank]['ignored']+=1
						else:
							result.append([1,ranks.getRankID(rank), valid_taxid, float(ab)])
							count[rank]['total']+=1
					else:
						count[rank]['ignored']+=1

		return np.array(result), count
		
	def parse_binning(input_file):
		try: # PANDAS implementation - faster
			import pandas as pd
			count = {'total':0,'ignored':0}
			header_count=0
			for line in open(input_file,'r'):
				if line[0]=="@" or line[0]=="#" or line[0]=="\n": header_count+=1
				else: break
			pandas_parsed = pd.read_csv(input_file, sep="\t", header=None, skiprows=header_count, names=('taxid','len'), usecols=[1,2], converters={'taxid': lambda txid: retrieveValidTaxID(int(txid))}, dtype={'len':int})
			count['total'] = pandas_parsed.shape[0]
			result = pandas_parsed[pandas_parsed.taxid!=0]
			count['ignored'] = count['total'] - result.shape[0]
			if verbose: 
				for row in pandas_parsed[pandas_parsed.taxid==0].itertuples(): print(("Ignored entry [%s] - Taxid not found" % (row.taxid)))
		except:	
			if verbose: print("(WARNING) pandas parsing failed")
			result = []
			count = {'total':0,'ignored':0}
			with open(input_file,'r') as f:
				for line in f:
					if line[0]=="@" or line[0]=="#" or line[0]=="\n": continue
					fields = line.rstrip().split('\t')
					taxid = int(fields[1])
					len = fields[2]
					count['total']+=1
					valid_taxid = retrieveValidTaxID(taxid)
					if valid_taxid: # Add taxid to results
						result.append([int(valid_taxid),int(len)])
					else:
						count['ignored']+=1
						if verbose: print(("Ignored entry [%s] - Taxid not found" % (taxid)))
		# result = [taxid, len]
		return np.array(result), count

	def b2p(binning_result):# Binning to Profiling
		# Sum read lengths for each entry (by taxid)
		leaf_len = defaultdict(int)
		for br in binning_result: leaf_len[br[0]]+=br[1]
	
		# From the counts (only leafs) reconstruct the taxonomic tree, summing up lengths to higher taxonomic levels
		for taxid,len in list(leaf_len.items()):
			#Start from the parent node because the entry itself (leaf) already contains the count
			txid = nodes[taxid]['parent']
			nlen = int(len)
			while txid!=1: # Until root node
				nlen = leaf_len[txid] + nlen # New length account for node accumulated length plus parent length
				leaf_len[txid] = nlen # Add new lenght to the parent
				txid = nodes[txid]['parent'] # Set next node
		
		# Generate the profile
		count = defaultdict(int)
		result = []
		for taxid,len in leaf_len.items():
			rank = nodes[taxid]['rank']
			# Only keep choosen ranks - They are ignored here but their counts were already summed up to the valid ranks
			if rank in ranks.ranks:
				count[rank]+=1
				result.append([1,ranks.getRankID(rank),taxid,float(len)])
		
		return np.array(result), count
	
	def isbioboxes(input_file):
		with open(input_file, 'r') as f: first_line = f.readline()
		return True if first_line.startswith("@") or first_line.startswith("#") else False
	
	#######################################################################
	
	if method=='db' or method=='p':
		isbiob = isbioboxes(input_file)
		if isbiob: print(" - %s (BioBoxes)" % input_file)
		else: print(" - %s (tsv)" % input_file)
		
		parsed_profile, profile_count = parse_profiling(input_file, isbiob)
		for rank in ranks.ranks: 
			print(("\t%s - %d entries (%d ignored)") % (rank, profile_count[rank]['total'], profile_count[rank]['ignored']))
			# Just warn user because lack of a rank can happen (when lineage is not well described)
			if profile_count[rank]['total']-profile_count[rank]['ignored']==0: print(("\t(WARNING) no valid entries found [%s]") % (rank))

	elif method=='b':
		if isbioboxes(input_file): print(" - %s (BioBoxes)" % input_file)
		else: print(" - %s (tsv)" % input_file)
		
		binning_result, binning_count = parse_binning(input_file)
		print(("\t%d lines (%d ignored)") % (binning_count['total'], binning_count['ignored']))
		
		parsed_profile, profile_count = b2p(binning_result)
		for rank in ranks.ranks: 
			print(("\t%s - %d entries") % (rank, profile_count[rank]))
			if profile_count[rank]==0: print(("\t(WARNING) no valid entries found [%s]") % (rank))

	#parsed_profile = np.array([Presence,RankID,TaxID,Val])
	return parsed_profile
