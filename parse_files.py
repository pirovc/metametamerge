import numpy as np
from collections import defaultdict

def parse_files(input_file, method, all_names_scientific, all_names_other, nodes, merged, ranks):
	
	def nameSearch(rank, name):
		## TODO - add name and create fake 'taxid'
		if (name,rank) in all_names_scientific:
			if len(all_names_scientific[(name,rank)])>1: print("(WARNING) ambiguous scientific name", rank, name)
			#all_names_scientific[(name,rank)][0] - first taxid (can be more than one if repeated entry)
			return all_names_scientific[(name,rank)][0]
		elif (name,rank) in all_names_other:
			if len(all_names_other[(name,rank)])>1: print("(WARNING) ambiguous name", rank, name)
			#all_names_other[(name,rank)][0] - first taxid (can be more than one if repeated entry)
			return all_names_other[(name,rank)][0]	
		else:	
			return 0

	def isbioboxes(input_file):
		with open(input_file, 'r') as f: first_line = f.readline()
		return True if first_line.startswith("@") or first_line.startswith("#") else False
	
	def parse_profiling_bioboxes(input_file):
		result = []
		with open(input_file,'r') as f:
			for line in f:
				if line[0]=="@" or line[0]=="#" or line[0]=="\n": continue
				fields = line.rstrip().split('\t')
				taxid = fields[0]
				rank = fields[1]
				name = fields[3].split("|")[-1]
				ab = fields[4]
				if rank not in ranks.ranks: #Verify rank only on profiling (binning is kept to account for the abundance estimation)
					print(("parse_files Rank not expected - ignoring entry [%s] rank [%s]" % (name,rank)))
					continue
				elif not taxid.isdigit(): #If taxid is not digit (ex: strain on CAMI data "103731.1")
					taxid = nameSearch(rank,name) #Try to find entry by name
				elif int(taxid) not in nodes: #If taxid is not present on nodes.dmp (ex: different version?)
					taxid = nameSearch(rank,name) #Try to find entry by name
				else: 
					taxid = int(taxid) # Valid taxid
				if taxid: result.append([1,ranks.getRankID(rank), taxid, float(ab)])
		return result
	
	def parse_profiling_tsv(input_file):
		result = []
		count = defaultdict(lambda: {'total':0,'ignored':0})
		with open(input_file,'r') as f:
			for line in f:
				taxid = 0
				rank, name_taxid, ab = line.rstrip().split('\t')
				count[rank]['total']+=1
				
				# Verify if rank is valid -> only on profiling (binning is kept to account for the abundance estimation)
				if rank not in ranks.ranks: 
					print(("Ignored entry [%s] rank [%s] - Rank not expected" % (name_taxid,rank)))
					count[rank]['ignored']+=1
					continue
					
				if name_taxid.isdigit(): #If taxid is not digit (ex: entry by name or strain on CAMI data "103731.1")
					name_taxid = int(name_taxid)
					if name_taxid in nodes: # Valid taxid
						taxid = name_taxid
					elif name_taxid in merged: # Search for updated taxid on merged.dmp
						print(("(WARNING) merged taxid [%d] rank [%s] -> [%s] rank [%s]") % (name_taxid, rank, merged[name_taxid], nodes[merged[name_taxid]]['rank']))
						# If found on merged.dmp, verify if rank is valid
						if nodes[merged[name_taxid]]['rank'] not in ranks.ranks: 
							print(("Ignored entry [%s] rank [%s] - Rank not expected" % (name_taxid,rank)))
							count[rank]['ignored']+=1
							continue
						else:
							taxid = merged[name_taxid]
					else:
						taxid = 0 # Not found
						
				# Try to find entry by name if taxid was not assigned directly (nodes.dmp, merged.dmp)	
				if not taxid: 
					taxid = nameSearch(rank,name_taxid) 
					if not taxid:
						print(("Ignored entry [%s] rank [%s] - Name not found" % (name_taxid,rank)))
						count[rank]['ignored']+=1
						continue
						
				result.append([1,ranks.getRankID(rank), taxid, float(ab)])
		
		print(count)
		return result, count
		
	def parse_binning_bioboxes(input_file):
		result = []
		with open(input_file,'r') as f:
			for line in f:
				if line[0]=="@" or line[0]=="#" or line[0]=="\n": continue
				fields = line.rstrip().split('\t')
				taxid = fields[1]
				len = fields[2]
				# Check if taxid is present on nodes.dmp
				if int(taxid) not in nodes:
					print(("parse_files Taxid not found - ignoring entry [%s] taxid [%s]" % (fields[0],taxid)))
				else:
					result.append([int(taxid),int(len)])
		return result

	def parse_binning_tsv(input_file):
		result = []
		with open(input_file,'r') as f:
			for line in f:
				fields = line.rstrip().split('\t')
				taxid = fields[1]
				len = fields[2]
				# Check if taxid is present on nodes.dmp
				if int(taxid) not in nodes:
					print(("parse_files Taxid not found - ignoring entry [%s] taxid [%s]" % (fields[0],taxid)))
				else:
					result.append([int(taxid),int(len)])
		return result
	
	def b2p_count(result,nodes):
		# Binning to Profiling
		
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
			if rank not in ranks.ranks:	
				pass#print("parse_files Rank not expected - ignoring entry on output [%s] rank [%s]" % (taxid,rank))
			else:
				r.append([1,ranks.getRankID(rank),taxid,float(len)])
		
		return r
	

	print(input_file)
	if method=='db' or method=='p':
		# Profiling parse_file rules:
		if isbioboxes(input_file):
			print("(BioBoxes)")
			r = parse_profiling_bioboxes(input_file)
		else:
			print("(tsv)")
			r, count = parse_profiling_tsv(input_file)
	elif method=='b':
		if isbioboxes(input_file):
			print("(BioBoxes)")
			r = b2p_count(parse_binning_bioboxes(input_file),nodes)
		else:
			print("(tsv)")
			r = b2p_count(parse_binning_tsv(input_file),nodes)
	
	for rank in ranks.ranks: print(rank, count[rank])
	#r = [Presence,RankID,TaxID,Val]
	return np.array(r)
