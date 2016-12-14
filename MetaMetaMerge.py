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

from metametamerge.Tools import Tools
from metametamerge.Databases import Databases
from metametamerge.Profile import Profile
from metametamerge.Ranks import Ranks

from metametamerge.parse_tax import parse_tax
from metametamerge.parse_files import parse_files

def hmean(l): return len(l) / sum(1. / val for val in l)

def main():
	version = '1.0'
	
	parser = argparse.ArgumentParser(description='MetaMetaMerge')
	parser.add_argument('-i', metavar='<input_files>', dest="input_files", nargs="*", help="Input binning or profiling files. Bioboxes or tsv format (see README)")
	parser.add_argument('-d', metavar='<database_profiles>', dest="database_profiles", nargs="*", help="Database profiles on the same order of the input files (see README)")
	parser.add_argument('-t', metavar='<tool_identifier>', dest="tool_identifier", type=str, help="Identifiers on the same order of the input files")
	parser.add_argument('-c', metavar='<tool_method>', dest="tool_method", type=str, help="Tools' method on the same order of the input files. p -> profiling / b -> binning")
	parser.add_argument('-n', metavar='<names_file>', dest="names_file", type=str, help="names.dmp from the NCBI Taxonomy database")
	parser.add_argument('-e', metavar='<nodes_file>', dest="nodes_file", type=str, help="nodes.dmp from the NCBI Taxonomy database")

	parser.add_argument('-b', metavar='<bins>', dest="bins", type=int, default=4,  help="Number of bins. Default: 4")
	parser.add_argument('-r', metavar='<cutoff>', dest="cutoff", type=float, default=0.0001, help="Minimum abundance/Maximum results for each taxonomic level (0: off / 0-1: minimum relative abundance / >=1: maximum number of identifications). Default: 0.0001") 	
	parser.add_argument('-f', metavar='<mode>', dest="mode", type=str, default="linear",  help="Result mode (precise, very-precise, linear, sensitive, very-sensitive). Default: linear")
	parser.add_argument('-s', metavar='<reversed_output>', dest="reversed_output", default=1, type=int, help="Consider only species level identifications and estimate upper taxonomic level identifications from it (0/1). Default: 1")
	
	parser.add_argument('-l', metavar='<detailed_output>', dest="detailed_output", default=0, type=int, help="Generate an additional detailed output (0/1). Default: 0")
	parser.add_argument('-o', metavar='<output_file>', dest="output_file", type=str, help="Output file")
	
	parser.add_argument('-v', action='version', version='%(prog)s ' + version)
	
	#parser.add_argument('-g', metavar='<ground_truth>', type=str, dest="ground_truth")
	args = parser.parse_args()

	ranks = Ranks(['superkingdom','phylum','class','order','family','genus','species'])
	output_folder = os.path.dirname(args.output_file)

	print("- - - - - - - - - - - - - - - - - - - - -")
	print("|\t\t\tMetaMetaMerge %s\t\t\t|" % version)
	print("- - - - - - - - - - - - - - - - - - - - -")
	print("Ranks: %s" % ranks.ranks)
	print("Input files: %s" % args.input_files)
	print("Database profiles: %s" % args.database_profiles)
	print("Identifiers: %s" % args.tool_identifier)
	print("Methods: %s" % args.tool_method)
	print("Names.dmp: %s" % args.names_file)
	print("Nodes.dmp: %s" % args.nodes_file)
	print("Mode: %s" % args.mode)
	print("Cutoff: %s" % args.cutoff)
	print("Bins: %s" % args.bins)
	print("Reversed: %s" % args.reversed_output)
	print("Detailed: %s" % args.detailed_output)
	print("Output file: %s" % args.output_file)
	print("- - - - - - - - - - - - - - - - - - - - -")

	# all_names_scientific, all_names_other -> defaultdict((name,rank): taxid})
	# nodes -> {taxid:{'parent':taxid,'rank':rank}} **** all nodes.dmp + names.dmp
	all_names_scientific, all_names_other, nodes = parse_tax(args.names_file, args.nodes_file, ranks)

	# Database profiles
	D = []
	dbs_count = defaultdict(int)
	for database_file in args.database_profiles:
		db = Databases(database_file, parse_files(database_file, 'db', all_names_scientific, all_names_other, nodes, ranks), ranks)
		D.append(db)
		# dbs_count -> {taxid: count}
		for taxid in db.getCol('TaxID'): dbs_count[taxid]+=1

	# Tools results
	T = []
	identifiers = args.tool_identifier.split(",")
	methods = args.tool_method.split(",")
	for idx,input_file in enumerate(args.input_files):
		T.append(Tools(input_file, identifiers[idx], methods[idx], parse_files(input_file, methods[idx], all_names_scientific, all_names_other, nodes, ranks), ranks, D[idx]))

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
	if args.cutoff>1:
		for tool in T: tool.filterMaxResults(int(args.cutoff))
	elif args.cutoff>0:
		for tool in T: tool.filterMinRelativeAbundance(args.cutoff)

	# Repeats due to same name,rank TODO - solve ambiguity
	for tool in T: tool.filterRepeatedNameRank()

	# Merged results
	merged = defaultdict(lambda: {'pres':0,'wepres':0,'ab':[]})
	tool_taxids = np.unique([taxid for tool in T for taxid in tool.getCol('TaxID')])
	for tool in T:
		tool.sort([('Abundance',-1)])
		for rankid,profilerank in tool:
			for pr in profilerank:
				# Only if has real abundance
				if pr['Abundance']==0:
					print("metametamerge WARNING zero abundance ", ranks.getRankName(rankid), nodes[taxid]['name'], tool.ident)
				else:
					taxid = pr['TaxID']
					if taxid in tool_taxids:
						# Presence
						merged[(taxid,rankid)]['pres'] += 1

						# # Weighted presence
						if taxid not in dbs_count:
							print("metametamerge ERROR entry not found in dbs_count ", ranks.getRankName(rankid), nodes[taxid]['name'], tool.ident)
							merged[(taxid,rankid)]['wepres'] = 0.5
						else:
							if merged[(taxid,rankid)]['pres'] > dbs_count[taxid]:
								merged[(taxid,rankid)]['wepres'] = ((dbs_count[taxid]+1)**2)/float(dbs_count[taxid])
								print("metametamerge ERROR more ident. than db elemnts ", nodes[taxid]['name'])
							else:
								merged[(taxid,rankid)]['wepres'] = ((merged[(taxid,rankid)]['pres']+1)**2)/float(dbs_count[taxid])

						# Abundances
						merged[(taxid,rankid)]['ab'].append(pr['Abundance'])
						
	pres = [val['wepres'] for (taxid,rankid),val in list(merged.items())]
	profile_merged = []
	bin_size = args.bins
	_, bin_edges = np.histogram(pres, bins=bin_size, range=(0,np.max(pres)))

	for (taxid,rankid),val in list(merged.items()):
		profile_merged.append([np.digitize([val['wepres']],bin_edges,right=True)[0],rankid,taxid,hmean(val['ab'])])

	profile_merged = Profile(np.array(profile_merged),ranks)

	#Sort by desc. presence and desc. abundance
	profile_merged.sort([('Abundance',-1),('Presence',-1)])

	# Apply cutoff based on the k presence (mode)
	profile_merged_kpres = []
	mode = {}
	max_pres = np.max(profile_merged.getCol('Presence'))
	for rankid,profilerank in profile_merged:
		for p in sorted(np.unique(profile_merged.getCol('Presence'))):
			# Organisms in the rank and n. of presence
			subset_pres = profilerank.getSubSet(profilerank.getCol('Presence')==p)
			x = subset_pres.getSize()

			if args.mode=="linear": #lin
				mode[(rankid, p)] = (x, p / float(max_pres))
			elif args.mode=="very-sensitive": #log_max
				mode[(rankid, p)] = (x, np.log(p + 3) / float(np.log(max_pres + 3)))
			elif args.mode=="sensitive": #log
				mode[(rankid, p)] = (x, np.log(p + 1) / float(np.log(max_pres + 1)))
			elif args.mode=="very-precise": #exp_max
				mode[(rankid, p)] = (x, (4 ** p) / float(4 ** max_pres))
			elif args.mode=="precise": #exp
				mode[(rankid, p)] = (x, (2 ** p) / float(2 ** max_pres))

			#mode[(rankid,p)] = (x,1) # NO CUTOFF
			#mode[(rankid,p)] = (x,1) if p==3 else (x,0)
			print("Rank: %s \t Presence: %s \t N.Orgs: %d \t Cutoff: %d (%.2f%%) " % (ranks.getRankName(rankid), p, x, mode[(rankid,p)][0]*mode[(rankid,p)][1], mode[(rankid,p)][1]*100))
			# For each subset selected
			for pr in subset_pres.getRow(0,mode[(rankid,p)][0]*mode[(rankid,p)][1]):
				profile_merged_kpres.append([pr['Presence'],rankid,pr['TaxID'],pr['Abundance']])
		print()

	# Add as a tool (+ normalize the abundance)
	profile_merged_kpres = Tools("", "merged", "p", np.array(profile_merged_kpres), ranks, "")

	# Sort merged results (ascending, based on position)
	profile_merged_kpres.sort([('Abundance',-1)])

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
		
	if args.reversed_output:
		### GENERATE LOWER RANKS FROM SPECIES
		reversed_ranks = defaultdict(lambda: {'Abundance':0,'Rank':''})
		for pr in profile_merged_kpres.profilerank[ranks.getRankID("species")]:
			txid = pr['TaxID']
			while txid!=1:
				if nodes[txid]['rank'] in ranks.ranks:
					reversed_ranks[txid]['Abundance']+=pr['Abundance']
					reversed_ranks[txid]['Rank']=nodes[txid]['rank']
				txid = nodes[txid]['parent']
		
		for r in ranks.ranks: #workaround for output sorted
			for txid, rp in list(reversed_ranks.items()):
				if rp['Rank']==r:
					out.write("%s\t%s\t%.16f\n" % (rp['Rank'],nodes[txid]['name'],rp['Abundance']))	
					if args.detailed_output: print_details(out_detailed,txid)
						
	else:
		for rankid, profilerank in profile_merged_kpres:
			for pr in profilerank:
				out.write("%s\t%s\t%.16f\n" % (ranks.getRankName(rankid),nodes[pr['TaxID']]['name'],pr['Abundance']))
				if args.detailed_output: print_details(out_detailed,pr['TaxID'])
				
	out.close()
	if args.detailed_output: out_detailed.close()
		
	#### PLOTS
	# if args.ground_truth:
		# import plots_dev
		# plots_dev.plots(T,all_names_scientific, all_names_other, nodes, profile_merged, profile_merged_kpres,mode,args.ground_truth,ranks,args.output_file)
    #### PLOTS
if __name__ == "__main__":
	main()
