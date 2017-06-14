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
import argparse, os
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

def hmean(l): return len(l) / sum(1. / val for val in l)

def main():
	version = '1.1'
	
	parser = argparse.ArgumentParser(description='MetaMetaMerge')
	parser.add_argument('-i', metavar='<input_files>', dest="input_files", nargs="*", help="Input binning or profiling files. Bioboxes or tsv format (see README)")
	parser.add_argument('-d', metavar='<database_profiles>', dest="database_profiles", nargs="*", help="Database profiles on the same order of the input files (see README)")
	parser.add_argument('-t', metavar='<tool_identifier>', dest="tool_identifier", type=str, help="Identifiers on the same order of the input files")
	parser.add_argument('-c', metavar='<tool_method>', dest="tool_method", type=str, help="Tools' method on the same order of the input files. p -> profiling / b -> binning")
	
	parser.add_argument('-n', metavar='<names_file>', dest="names_file", type=str, help="names.dmp from the NCBI Taxonomy database")
	parser.add_argument('-e', metavar='<nodes_file>', dest="nodes_file", type=str, help="nodes.dmp from the NCBI Taxonomy database")
	parser.add_argument('-m', metavar='<merged_file>', dest="merged_file", type=str, help="merged.dmp from the NCBI Taxonomy database")

	parser.add_argument('-b', metavar='<bins>', dest="bins", type=int, default=4,  help="Number of bins. Default: 4")
	parser.add_argument('-r', metavar='<cutoff>', dest="cutoff", type=float, default=0.0001, help="Minimum abundance/Maximum results for each taxonomic level (0: off / 0-1: minimum relative abundance / >=1: maximum number of identifications). Default: 0.0001") 	
	parser.add_argument('-f', metavar='<mode>', dest="mode", type=str, default="linear",  help="Result mode (precise, very-precise, linear, sensitive, very-sensitive, no-cutoff). Default: linear")
	parser.add_argument('-s', metavar='<ranks>', dest="ranks", default="species", type=str, help="Consider only the specified taxonomic ranks (comma separated) and estimate upper taxonomic levels from it. Supressing this  option will merge every rank independently for every tool. (superkingdom, phylum , class, order, family, genus, species, off). Default: species")
	
	parser.add_argument('--verbose', action='store_true')
	parser.add_argument('-l', metavar='<detailed_output>', dest="detailed_output", default=0, type=int, help="Generate an additional detailed output (0/1). Default: 0")
	parser.add_argument('-o', metavar='<output_file>', dest="output_file", type=str, help="Output file")
	
	parser.add_argument('-v', action='version', version='%(prog)s ' + version)
	
	args = parser.parse_args()

	all_ranks = Ranks(['superkingdom','phylum','class','order','family','genus','species'])
	if args.ranks:
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
	print("Mode: %s" % args.mode)
	print("Cutoff: %s" % args.cutoff)
	print("Bins: %s" % args.bins)
	print("Ranks: %s" % ', '.join(ranks.ranks))
	print("Verbose: %s" % args.verbose)
	print("Detailed: %s" % args.detailed_output)
	print("Output file: %s" % args.output_file)
	print("- - - - - - - - - - - - - - - - - - - - -")
	print()

	print("Parsing taxonomy (names, nodes, merged) ... ")
	# all_names_scientific, all_names_other -> defaultdict((name,rank): taxid})
	# nodes -> {taxid:{'parent':taxid,'rank':rank}} **** all nodes.dmp + names.dmp
	all_names_scientific, all_names_other, nodes, merged = parse_tax(args.names_file, args.nodes_file, args.merged_file, ranks)
	print()
	
	print("Reading database profiles ...")
	# Database profiles
	D = []
	dbs_count = defaultdict(int)
	for database_file in args.database_profiles:
		db = Databases(database_file, parse_files(database_file, 'db', all_names_scientific, all_names_other, nodes, merged, ranks, args.verbose), ranks)
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
		tool = Tools(input_file, identifiers[idx], methods[idx], parse_files(input_file, methods[idx], all_names_scientific, all_names_other, nodes, merged, ranks, args.verbose), ranks, args.verbose, D[idx])
		T.append(tool)
	print()
	
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
	if args.cutoff>=1:
		print("Filtering relative abundances (= 0)")
		for tool in T: tool.filterMinRelativeAbundance(0, ranks)
		print()
		print("Filtering profiles (max. results = %d) ..." % args.cutoff)
		for tool in T: tool.filterMaxResults(int(args.cutoff), ranks)
		print()
	elif args.cutoff>0:
		print("Filtering profiles (min. relative abundance <= %f) ..." % args.cutoff)
		for tool in T: tool.filterMinRelativeAbundance(args.cutoff, ranks)
		print()
	else:
		print("Filtering relative abundances (= 0)")
		for tool in T: tool.filterMinRelativeAbundance(0, ranks)
		print()
		
	# Merged results
	print("Merging profiles ...")
	merged = defaultdict(lambda: {'Presence':0,'Score':0,'Abundance':[]})
	for tool in T:
		for rankid,profilerank in tool:
			for pr in profilerank:
				taxid = pr['TaxID']
				merged[(taxid,rankid)]['Presence'] += 1 # Counting pr['TaxID'] would be duplicated for merged
				merged[(taxid,rankid)]['Score'] = ((merged[(taxid,rankid)]['Presence']+1)**2)/float(dbs_count[taxid]+1)
				merged[(taxid,rankid)]['Abundance'].append(pr['Abundance'])
	
	# List of scores
	scores = [val['Score'] for (taxid,rankid),val in list(merged.items())]
	# Divide scores range into bin groups 
	_, bin_edges = np.histogram(scores, bins=args.bins, range=(0,len(T)+1))
	
	# Create an Profile, setting the bin number (1...args.bins) according to the edges and taking harmonic mean out of the abundances
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

			print("\tbin: %d \t # taxons: %d (%d kept - %.2f%%) " % (bin, bin_c, mode[(rankid,bin)][0]*mode[(rankid,bin)][1], mode[(rankid,bin)][1]*100))
				
			if bin_c:
				# Add entry only for taxons above the cutoff
				for pr in profile_bin.getRow(0,int(mode[(rankid,bin)][0]*mode[(rankid,bin)][1])):
					profile_merged_mode.append([pr['Presence'],rankid,pr['TaxID'],pr['Abundance']])
				
	# Add as a tool (+ normalize the abundance)
	profile_merged_mode = Tools("", "merged", "p", np.array(profile_merged_mode), ranks, args.verbose, "")

	# Sort merged results (ascending, based on position)
	profile_merged_mode.sort([('Abundance',-1)])
	
	print()
	print("Final merged profile:")
	for rankid,profilerank in profile_merged_mode:
		print(("\t%s - %d entries") % (ranks.getRankName(rankid),profilerank.getSize()))
	print()
	

	
	# Print merged results
	out = open(args.output_file,'w')
	if args.detailed_output:
		out_detailed = open(args.output_file + ".detailed",'w')
		out_detailed.write("rank\ttaxid\tname\tmetametamerge\t")
		out_detailed.write('\t'.join([t.ident for t in T]))
		out_detailed.write('\n')
		
	def print_details(od,tid):
		od.write("%s\t%d\t%s\t%.16f\t" % (rp['Rank'],tid,nodes[tid]['name'],rp['Abundance']))
		for toolid,tool in enumerate(T):
			pres_t = tool.getSubSet('TaxID',tid).getCol('Abundance')
			if pres_t: #Identified by the tool
				od.write("%.16f\t" % np.sum(pres_t))
			else:
				pres_d = D[toolid].getSubSet('TaxID',tid).getCol('Abundance')
				if pres_d: # Not identified but present 
					od.write("0\t")
				else:# Not present in the db
					od.write("-1\t")
		od.write("\n")
		
	if args.ranks:
		### GENERATE LOWER RANKS FROM SPECIES
		reversed_ranks = defaultdict(lambda: {'Abundance':0,'Rank':''})
		for pr in profile_merged_mode.profilerank[ranks.getRankID("species")]:
			txid = pr['TaxID']
			while txid!=1:
				if nodes[txid]['rank'] in all_ranks.ranks:
					reversed_ranks[txid]['Abundance']+=pr['Abundance']
					reversed_ranks[txid]['Rank']=nodes[txid]['rank']
				txid = nodes[txid]['parent']
		
		for r in all_ranks.ranks: #workaround for output sorted
			for txid, rp in list(reversed_ranks.items()):
				if rp['Rank']==r:
					out.write("%s\t%s\t%.16f\n" % (rp['Rank'],nodes[txid]['name'],rp['Abundance']))	
					if args.detailed_output: print_details(out_detailed,txid)
						
	else:
		for rankid, profilerank in profile_merged_mode:
			for pr in profilerank:
				out.write("%s\t%s\t%.16f\n" % (ranks.getRankName(rankid),nodes[pr['TaxID']]['name'],pr['Abundance']))
				if args.detailed_output: print_details(out_detailed,pr['TaxID'])
				
	out.close()
	if args.detailed_output: out_detailed.close()
		
	#### PLOTS
	# if args.ground_truth:
		# import plots_dev
		# plots_dev.plots(T,all_names_scientific, all_names_other, nodes, profile_merged, profile_merged_mode,mode,args.ground_truth,ranks,args.output_file)
    #### PLOTS
if __name__ == "__main__":
	main()
