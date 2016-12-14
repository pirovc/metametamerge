from Tools import Tools
from Databases import Databases
from Profile import Profile
from Ranks import Ranks
from parse_tax import parse_tax
from parse_files import parse_files

from collections import defaultdict, OrderedDict
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def plots(T,all_names_scientific, all_names_other ,  nodes, profile_merged,profile_merged_kpres,cutoff,ground_truth,ranks,out_file):
	rank = "species"	
	G = Tools(ground_truth, 'gt', 'p', parse_files(ground_truth, 'p', all_names_scientific, all_names_other, nodes, ranks), ranks)
	gt_profile = Profile([],[],[G.profilerank[ranks.getRankID(rank)]])

	# Plots
	fig, ax = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False)

	# 1
	adj = 0
	old_pres = 0
	tp = defaultdict(int)
	fp = defaultdict(int)
	for rankid, profilerank in profile_merged:
		# Just one rank - Compare string (from manual parameter)
		if rankid==ranks.getRankID(rank):
			for p,pr in enumerate(profilerank):
				if old_pres==0: old_pres=pr['Presence']
				if pr['TaxID'] in gt_profile.getCol('TaxID'):
					c = 'g' 
					tp[pr['Presence']]+=1
				else:
					c = 'r' 
					fp[pr['Presence']]+=1
				if pr['Presence']!=old_pres: 
					co = cutoff[(rankid,old_pres)][0]*cutoff[(rankid,old_pres)][1]
					ax[0,0].plot([old_pres-0.5,old_pres+0.5],[co,co],lw=2,color="black")	
					adj = p-1
				ax[0,0].plot([pr['Presence']-0.4,pr['Presence']+0.4],[p-adj,p-adj],color=c,lw=1)
				old_pres=pr['Presence']
				#print p, pr['Presence'], p-adj
			co = cutoff[(rankid,old_pres)][0]*cutoff[(rankid,old_pres)][1]
			ax[0,0].plot([old_pres-0.5,old_pres+0.5],[co,co],lw=2,color="black")	
	ax[0,0].set_xlabel('Presence [' + rank + ']')
	ax[0,0].set_ylabel('Ordered output')

	#################
	tool_short = []
	#sortT = iter(sorted(T.profile.items()))
	for id,tool in enumerate(T):
		tool_short.append(tool.ident)
		tool.sort([('Abundance',-1)])
		for rankid, profilerank in tool:
			# Just one rank - Compare string (from manual parameter)
			if rankid==ranks.getRankID(rank):
				for p,pr in enumerate(profilerank):
					# Compare string because uniqueids where generated differently
					if pr['TaxID'] in gt_profile.getCol('TaxID'):
						c = 'g'
					else:
						c = 'r'
					ax[0,1].plot([id-0.4,id+0.4],[p,p],color=c,lw=1)

	tool_short.append('MERGED')
	for rankid, profilerank in profile_merged_kpres:
		# Just one rank - Compare string (from manual parameter)
		if rankid==ranks.getRankID(rank):
			for p,pr in enumerate(profilerank):
				# Compare string because uniqueids where generated differently
				if pr['TaxID'] in gt_profile.getCol('TaxID'):
					c = 'g'
				else:
					c = 'r'
				ax[0,1].plot([id+1-0.4,id+1+0.4],[p,p],color=c,lw=1)
	ax[0,1].set_xticks([r for r in range(len(tool_short))])
	ax[0,1].set_xticklabels(tool_short, rotation=90)
	ax[0,1].set_ylabel('Ordered output')

	#################
	def plotTool(prof,x,y,x_lim=None,y_lim=None):
		pres_tool = {tool.file:{pres:{'tp':0,'fp':0} for pres in np.unique(profile_merged.getCol('Presence'))} for tool in T}
		for rankid, profilerank in prof:
			# Just one rank - Compare string (from manual parameter)
			if rankid==ranks.getRankID(rank):
				for p,pr in enumerate(profilerank):
					for tool in T:
						if pr['TaxID'] in tool.profilerank[rankid].getCol('TaxID'):
							if pr['TaxID'] in gt_profile.getCol('TaxID'):
								pres_tool[tool.file][pr['Presence']]['tp']+=1
							else:
								pres_tool[tool.file][pr['Presence']]['fp']+=1	
		#print pres_tool
		colors = ['k','g','r','c','m','y','b']
		for p in sorted(tp.keys()):
			for id,(tool,pres) in enumerate(sorted(pres_tool.items()),1):
				ax[x,y].bar(p + (id/float(len(pres_tool)+1)),pres[p]['tp'],1/float(len(pres_tool)+1), color=colors[id], label=tool)
				ax[x,y].bar(p + (id/float(len(pres_tool)+1)),-pres[p]['fp'],1/float(len(pres_tool)+1), color=colors[id], label=tool)
		ax[x,y].axhline(0, color='black')
		
		ax[x,y].set_ylabel('FP | TP')
		handles, labels = ax[x,y].get_legend_handles_labels()
		if x_lim: ax[x,y].set_xlim(x_lim)
		if y_lim: ax[x,y].set_ylim(y_lim)
		#ax[x,y].legend([handles[r] for r in range(0,len(pres_tool)*2,2)],[labels[r].split("/")[-1].split("_")[0] for r in range(0,len(pres_tool)*2,2)],bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=len(pres_tool)/2, mode="expand", borderaxespad=0.)

		return ax[x,y].get_xlim(), ax[x,y].get_ylim()
	x_lim,y_lim = plotTool(profile_merged,1,0)
	ax[1,0].set_xlabel('Presence [' + rank + '] before filter')
	plotTool(profile_merged_kpres,1,1,x_lim,y_lim)
	ax[1,1].set_xlabel('Presence [' + rank + '] after filter')

	
	#plt.show()
	fig.set_size_inches(15,8)
	plt.savefig(out_file + ".png",dpi=300)
