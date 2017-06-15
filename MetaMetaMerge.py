#!/usr/bin/python3
# The MIT License (MIT)
# 
# Copyright (c) 2016 - Vitor C. Piro - PiroV@rki.de - vitorpiro@gmail.com
# Robert Koch-Institut, Germany
# All rights reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import numpy as np
np.set_printoptions(suppress=True, precision=16, threshold=10000000)
import argparse, os, math
from collections import defaultdict

# from metametamerge.Tools import Tools
# from metametamerge.Databases import Databases
# from metametamerge.Profile import Profile
# from metametamerge.Ranks import Ranks
# from metametamerge.parse_tax import parse_tax
# from metametamerge.parse_files import parse_files
from Tools import Tools
from Databases import Databases
from Profile import Profile
from Ranks import Ranks
from parse_tax import parse_tax
from parse_files import parse_files

def main():
	version = '1.1'
	
	parser = argparse.ArgumentParser(description='MetaMetaMerge')
	parser.add_argument('-i', '--input-files', metavar='<input_files>', dest="input_files", nargs="*", required=True, help="Input (binning or profiling) files. Bioboxes or tsv format (see README)")
	parser.add_argument('-d', '--database-profiles', metavar='<database_profiles>', dest="database_profiles", nargs="*", required=True, help="Database profiles on the same order of the input files (see README)")
	parser.add_argument('-t', '--tool-identifier', metavar='<tool_identifier>', dest="tool_identifier", type=str, required=True, help="Comma-separated identifiers on the same order of the input files")
	parser.add_argument('-c', '--tool-method', metavar='<tool_method>', dest="tool_method", type=str, required=True, help="Comma-separated methods on the same order of the input files (p -> profiling / b -> binning)")
	
	parser.add_argument('-n', '--names-file', metavar='<names_file>', dest="names_file", type=str, required=True, help="names.dmp from the NCBI Taxonomy database")
	parser.add_argument('-e', '--nodes-file', metavar='<nodes_file>', dest="nodes_file", type=str, required=True, help="nodes.dmp from the NCBI Taxonomy database")
	parser.add_argument('-m', '--merged-file', metavar='<merged_file>', dest="merged_file", type=str, required=True, help="merged.dmp from the NCBI Taxonomy database")

	parser.add_argument('-b', '--bins', metavar='<bins>', dest="bins", type=int, default=4, help="Number of bins. Default: 4")
	parser.add_argument('-r', '--cutoff', metavar='<cutoff>', dest="cutoff", type=float, default=0.0001, help="Minimum abundance/Maximum results for each taxonomic level (0: off / 0-1: minimum relative abundance / >=1: maximum number of identifications). Default: 0.0001")
	parser.add_argument('-f', '--mode', metavar='<mode>', dest="mode", type=str, default="linear",  help="Result mode (precise, very-precise, linear, sensitive, very-sensitive, no-cutoff). Default: linear")
	parser.add_argument('-s', '--ranks', metavar='<ranks>', dest="ranks", default="species", type=str, help="Comma-separated list of ranks to be independently merged (superkingdom,phylum,class,order,family,genus,species,all). Default: species")
	
	parser.add_argument('-o', '--output-file', metavar='<output_file>', dest="output_file", type=str, required=True, help="Output file")
	parser.add_argument('-p', '--output-type', metavar='<output_type>', dest="output_type", default="bioboxes", type=str, help="Output type (tsv, bioboxes). Default: bioboxes")
	parser.add_argument('--detailed', action='store_true', dest="detailed", help="Generate an additional detailed output with individual normalized abundances for each tool, where: 0 -> not identified but present in the database, -1 not present in the database.")
	parser.add_argument('--verbose', action='store_true', dest="verbose", help="Verbose output log")
	
	parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version)
	
	args = parser.parse_args()

	all_ranks = Ranks(['superkingdom','phylum','class','order','family','genus','species'])
	if args.ranks and args.ranks!="all":
		ranks = Ranks([r.strip() for r in args.ranks.split(",")])
	else:
		ranks = all_ranks

	output_folder = os.path.dirname(args.output_file)

	print("- - - - - - - - - - - - - - - - - - - - -")
	print("           MetaMetaMerge %s" % version)
	print("- - - - - - - - - - - - - - - - - - - - -")
	print("Input files: ")
	for i,file in enumerate(args.input_files):
		print((" %s (%s) %s %s") % (args.tool_identifier.split(",")[i],args.tool_method.split(",")[i],file,args.database_profiles[i]))
	print("Taxonomy: \n %s, %s, %s" % (args.names_file,args.nodes_file,args.merged_file))
	print("Bins: %s" % args.bins)
	print("Cutoff: %s" % args.cutoff)
	print("Mode: %s" % args.mode)
	print("Ranks: %s" % ', '.join(ranks.ranks))
	print("Output file (type): %s (%s)" % (args.output_file,args.output_type))
	print("Verbose: %s" % args.verbose)
	print("Detailed: %s" % args.detailed)
	print("- - - - - - - - - - - - - - - - - - - - -")

	print()
	print("Parsing taxonomy (names, nodes, merged) ... ")
	# all_names_scientific, all_names_other -> defaultdict((name,rank): taxid})
	# nodes -> {taxid:{'parent':taxid,'rank':rank,'name':name}} **** all nodes.dmp + names.dmp
	all_names_scientific, all_names_other, nodes, merged = parse_tax(args.names_file, args.nodes_file, args.merged_file, all_ranks)
	
	print()
	print("Reading database profiles ...")
	# Database profiles
	D = []
	dbs_count = defaultdict(int)
	for database_file in args.database_profiles:
		db = Databases(database_file, parse_files(database_file, 'db', all_names_scientific, all_names_other, nodes, merged, ranks, args.verbose), ranks)		
		# Merge repeated taxids (when taxid changes in the new taxonomy version)
		merged_taxids = db.mergeRepeatedTaxIDs()
		if merged_taxids: print(("\t%d taxons with merged entries [%s]") % (len(merged_taxids),",".join([str(int(taxid)) for taxid in merged_taxids])))
		print(("\tTotal - %d taxons") % (db.getSize()))
		D.append(db)
		# dbs_count -> {taxid: count}
		for taxid in db.getCol('TaxID'): dbs_count[taxid]+=1
	
	print()
	print("Reading profiles ...")
	# Tools results
	T = []
	identifiers = args.tool_identifier.split(",")
	methods = args.tool_method.split(",")
	for idx,input_file in enumerate(args.input_files):
		tool = Tools(input_file, identifiers[idx], methods[idx], parse_files(input_file, methods[idx], all_names_scientific, all_names_other, nodes, merged, ranks, args.verbose), ranks, args.verbose)

		# Check presence on it's on database profile
		tool.checkDB(D[idx], ranks, args.verbose)
		
		# Merge repeated taxids (when taxid changes in the new taxonomy version)		
		merged_taxids = tool.mergeRepeatedTaxIDs()
		if merged_taxids: print(("\t%d taxons with merged entries [%s]") % (len(merged_taxids),",".join([str(int(taxid)) for taxid in merged_taxids])))
		
		# Estimate abundance for binning methods
		if methods[idx]=='b': tool.estimateAbundance(D[idx], args.verbose)
		
		print(("\tTotal - %d taxons") % (tool.getSize()))
		T.append(tool)

	# Print tool output for plots (before cutoff - only with normalized/estimated abundances - without entries not found in the DB!!)
	# if args.ground_truth:
		# for tool in T:
			# out = open(output_folder + "/metameta_" + tool.ident + ".out",'w')
			# tool.sort([('Abundance',-1)])
			# for rankid, profilerank in tool:
				# for pr in profilerank:
					# out.write("%s\t%s\t%.16f\n" % (ranks.getRankName(rankid),int(pr['TaxID']),pr['Abundance']))
			# out.close()

	# Filter max results
	print()
	if args.cutoff>=1:
		print("Filtering relative abundances (= 0)")
		for tool in T: tool.filterMinRelativeAbundance(0, ranks)
		print()
		print("Filtering profiles (max. results = %d) ..." % args.cutoff)
		for tool in T: tool.filterMaxResults(int(args.cutoff), ranks)
	elif args.cutoff>0:
		print("Filtering profiles (min. relative abundance <= %f) ..." % args.cutoff)
		for tool in T: tool.filterMinRelativeAbundance(args.cutoff, ranks)
	else:
		print("Filtering relative abundances (= 0)")
		for tool in T: tool.filterMinRelativeAbundance(0, ranks)
		
	# Merged results
	print()
	print("Merging profiles ...")
	merged = defaultdict(lambda: {'Presence':0,'Score':0,'Abundance':[]})
	for tool in T:
		for rankid,profilerank in tool:
			for pr in profilerank:
				taxid = pr['TaxID']
				merged[(taxid,rankid)]['Presence'] += 1 # Counting 1 instead of pr['TaxID'] since it could be >1 for merged entries
				merged[(taxid,rankid)]['Score'] = ((merged[(taxid,rankid)]['Presence']+1)**2)/float(dbs_count[taxid]+1)
				merged[(taxid,rankid)]['Abundance'].append(pr['Abundance'])
	
	# List of scores
	scores = [val['Score'] for (taxid,rankid),val in list(merged.items())]
	# Divide scores range into bin groups 
	_, bin_edges = np.histogram(scores, bins=args.bins, range=(0,len(T)+1))
	
	# Create an Profile, setting the bin number (1...args.bins) according to the edges and taking harmonic mean out of the abundances
	def hmean(l): return len(l) / sum(1. / val for val in l)
	profile_merged = []
	for (taxid,rankid),val in list(merged.items()):
		profile_merged.append([np.digitize([val['Score']],bin_edges,right=True)[0],rankid,taxid,hmean(val['Abundance'])])
	profile_merged = Profile(np.array(profile_merged),ranks)

	#Sort by desc. presence and desc. abundance
	profile_merged.sort([('Abundance',-1),('Presence',-1)])
	
	for rankid,profilerank in profile_merged:
		print(("\t%s - %d entries") % (ranks.getRankName(rankid),profilerank.getSize()))
	
	print()
	print("Aplying guided cutoff ...")
	# Apply guided cutoff based on pre-defined functions (mode)
	profile_merged_mode = []
	mode = {}
	bin_n = len(bin_edges)-1 #args.bins
	for rankid,profilerank in profile_merged:
		print("\t%s" % ranks.getRankName(rankid))
		for bin in range(1,bin_n+1):
			profile_bin = profilerank.getSubSet(profilerank.getCol('Presence')==bin)
			bin_c = profile_bin.getSize()
			
			if args.mode=="linear": #lin
				mode[(rankid, bin)] = (bin_c, bin / float(bin_n))
			elif args.mode=="very-sensitive": #log_max
				mode[(rankid, bin)] = (bin_c, np.log(bin + 3) / float(np.log(bin_n + 3)))
			elif args.mode=="sensitive": #log
				mode[(rankid, bin)] = (bin_c, np.log(bin + 1) / float(np.log(bin_n + 1)))
			elif args.mode=="very-precise": #exp_max
				mode[(rankid, bin)] = (bin_c, (4 ** bin) / float(4 ** bin_n))
			elif args.mode=="precise": #exp
				mode[(rankid, bin)] = (bin_c, (2 ** bin) / float(2 ** bin_n))
			elif args.mode=="no-cutoff": #do not apply cut-off, just merge results
				mode[(rankid, bin)] = (bin_c, 1) 

			print("\tbin: %d \t # taxons: %d (%d kept - %.2f%%) " % (bin, bin_c, math.ceil(mode[(rankid,bin)][0]*mode[(rankid,bin)][1]), mode[(rankid,bin)][1]*100))
				
			if bin_c:
				# Add entry only for taxons above the cutoff
				for pr in profile_bin.getRow(0,math.ceil(mode[(rankid,bin)][0]*mode[(rankid,bin)][1])):
					profile_merged_mode.append([pr['Presence'],rankid,pr['TaxID'],pr['Abundance']])
	
	# Add as a tool (+ normalize the abundance)
	profile_merged_mode = Tools("", "merged", "p", np.array(profile_merged_mode), ranks, args.verbose)
	
	# If only a subset of ranks were analyzed, estimate missing ranks
	if ranks.ranks!=all_ranks.ranks:
		# Connection between ranks and all_ranks should be done by name (rankid is different)	
		highest_rank_all_ranks = max([all_ranks.getRankID(rank_name) for rank_name in ranks.ranks if rank_name in all_ranks.ranks])
		missing_ranks = [rank_name for rank_name in all_ranks.ranks if rank_name not in ranks.ranks and all_ranks.getRankID(rank_name) < highest_rank_all_ranks]
		print()
		print(("Estimating missing ranks (%s) from %s ...") % (",".join(missing_ranks),all_ranks.getRankName(highest_rank_all_ranks)))
		
		reversed_ranks = defaultdict(lambda: {'Abundance':0,'Rank':''})
		for pr in profile_merged_mode.profilerank[ranks.getRankID(all_ranks.getRankName(highest_rank_all_ranks))]:
			txid = pr['TaxID']
			while txid!=1:
				if nodes[txid]['rank'] in missing_ranks:
					reversed_ranks[txid]['Abundance']+=pr['Abundance']
					reversed_ranks[txid]['Rank']=nodes[txid]['rank']
				txid = nodes[txid]['parent']
		
		profile_estimated = []
		# Add previous entries
		for rankid, profilerank in profile_merged_mode:
			for pr in profilerank:
				# Covert rankid (ranks) to rankid(all_ranks)
				profile_estimated.append([pr['Presence'],all_ranks.getRankID(ranks.getRankName(rankid)),pr['TaxID'],pr['Abundance']])
		# Add estimated entries
		for txid, rr in list(reversed_ranks.items()):
			profile_estimated.append([1,all_ranks.getRankID(rr['Rank']),txid,rr['Abundance']])
		
		# Add as a tool (+ normalize the abundance)
		profile_merged_mode = Tools("", "merged", "p", np.array(profile_estimated), all_ranks, args.verbose)
	
	# Sort merged results (ascending, based on position)
	profile_merged_mode.sort([('Abundance',-1)])
	print()
	print("Final merged profile:")
	for rankid,profilerank in profile_merged_mode:
		print(("\t%s - %d entries") % (all_ranks.getRankName(rankid),profilerank.getSize()))

	# Print final merged profile
	out = open(args.output_file,'w')
	if args.output_type=="tsv":
		for rankid, profilerank in profile_merged_mode:
			for pr in profilerank:
				out.write("%s\t%s\t%.16f\n" % (all_ranks.getRankName(rankid),nodes[pr['TaxID']]['name'],pr['Abundance']))	
	else: #bioboxes
		out.write("# Taxonomic Profiling Output\n")
		out.write("@SampleID:%s\n" % args.output_file)
		out.write("@Version:0.9.3\n")
		out.write("@Ranks:%s\n" % '|'.join(all_ranks.ranks))
		out.write("@TaxonomyID:%s\n" % args.nodes_file)
		out.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
		for rankid, profilerank in profile_merged_mode:
			for pr in profilerank:
				lineage = defaultdict(int)
				txid = int(pr['TaxID'])
				while txid!=1:
					lineage[nodes[txid]['rank']] = txid
					txid = nodes[txid]['parent']
				name_lineage = []
				taxid_lineage = []
				for r in all_ranks.ranks[0:rankid+1]:
					name_lineage.append(nodes[lineage[r]]['name'] if lineage[r] else "")
					taxid_lineage.append(str(lineage[r]) if lineage[r] else "")
				txid = int(pr['TaxID'])
				out.write(("%d\t%s\t%s\t%s\t%.6f\n" % (
						txid, 
						nodes[txid]['rank'],
						"|".join(taxid_lineage),
						"|".join(name_lineage),
						pr['Abundance']
						)))
	out.close()

	# Print detailed profile
	if args.detailed:
		out_detailed = open(args.output_file + ".detailed",'w')
		out_detailed.write("rank\ttaxid\tname\tmetametamerge\t")
		out_detailed.write('\t'.join([t.ident for t in T]))
		out_detailed.write('\n')
		for rankid, profilerank in profile_merged_mode:
			for pr in profilerank:
				out_detailed.write("%s\t%d\t%s\t%.16f\t" % (all_ranks.getRankName(rankid),pr['TaxID'],nodes[pr['TaxID']]['name'],pr['Abundance']))
				for toolid,tool in enumerate(T):
					pres_t = tool.getSubSet('TaxID',pr['TaxID']).getCol('Abundance')
					if pres_t: #Identified by the tool
						out_detailed.write("%.16f\t" % pres_t[0])
					else:
						pres_d = D[toolid].getSubSet('TaxID',pr['TaxID']).getCol('Abundance')
						if pres_d: # Not identified but present 
							out_detailed.write("0\t")
						else:# Not present in the db
							out_detailed.write("-1\t")
				out_detailed.write("\n")
		out_detailed.close()

if __name__ == "__main__":
	main()
