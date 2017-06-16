import sys
import argparse, math
import numpy as np
np.set_printoptions(suppress=True, threshold=10000000)

import matplotlib as mpl
mpl.use('Agg') #to plot without X

import matplotlib.pyplot as plt
from collections import defaultdict

from Tools import Tools
from Profile import Profile
from parse_files import parse_files
from parse_tax import parse_tax
from Ranks import Ranks

def main():	
	parser = argparse.ArgumentParser(description='metameta')
	parser.add_argument('-i', metavar='<input_files>', nargs="*", dest="input_files", help="...")
	parser.add_argument('-t', metavar='<tool_identifier>', dest="tool_identifier", help="...")
	parser.add_argument('-c', metavar='<tool_method>', dest="tool_method", help="...")
	parser.add_argument('-g', metavar='<ground_truth>', type=str, dest="ground_truth", help="...")
	parser.add_argument('-r', metavar='<rank>', type=str, dest="rank", default="species", help="...")
	parser.add_argument('-n', '--names-file', metavar='<names_file>', dest="names_file", type=str, required=True, help="names.dmp from the NCBI Taxonomy database")
	parser.add_argument('-e', '--nodes-file', metavar='<nodes_file>', dest="nodes_file", type=str, required=True, help="nodes.dmp from the NCBI Taxonomy database")
	parser.add_argument('-m', '--merged-file', metavar='<merged_file>', dest="merged_file", type=str, required=True, help="merged.dmp from the NCBI Taxonomy database")
	parser.add_argument('-o', metavar='<output_plot>', type=str, dest="output_plot", help="...")
	args = parser.parse_args()
	
	ranks = Ranks(['superkingdom','phylum','class','order','family','genus','species'])
	
	all_names_scientific, all_names_other, nodes, merged = parse_tax(args.names_file, args.nodes_file, args.merged_file, ranks)

	# Tools
	T = []
	identifiers = args.tool_identifier.split(",")
	methods = args.tool_method.split(",")
	for idx,input_file in enumerate(args.input_files):
		T.append(Tools(input_file, identifiers[idx], methods[idx], parse_files(input_file, methods[idx], all_names_scientific, all_names_other, nodes, merged, ranks, False), ranks, False))
		
	# Ground truth
	G = Tools(args.ground_truth, 'gt', 'p', parse_files(args.ground_truth, 'p', all_names_scientific, all_names_other, nodes, merged, ranks, False), ranks, False)
	G.normalizeAbundance()
	
	# Plots	
	fig, ax = plt.subplots(nrows=3, ncols=3, sharex=False, sharey=False)

	#('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd')
	markers = ('>','8','d','p','s','D','x','*')
	marker_size = 6
	colors = ("#000000","#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#FEFFE6","#1B4400","#4FC601","#3B5DFF","#4A3B53","#FF2F80","#61615A","#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA","#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299","#300018","#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C","#6F0062","#0CBD66","#EEC3FF","#456D75","#B77B68","#7A87A1","#788D66","#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459","#456648","#0086ED","#886F4C","#34362D","#B4A8BD","#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F","#938A81","#575329","#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757","#C8A1A1","#1E6E00","#7900D7","#A77500","#6367A9","#A05837","#6B002C","#772600","#D790FF","#9B9700","#549E79","#FFF69F","#201625","#72418F","#BC23FF","#99ADC0","#3A2465","#922329","#5B4534","#FDE8DC","#404E55","#0089A3","#CB7E98","#A4E804","#324E72","#6A3A4C","#83AB58","#001C1E","#D1F7CE","#004B28","#C8D0F6","#A3A489","#806C66","#222800","#BF5650","#E83000","#66796D","#DA007C","#FF1A59","#8ADBB4","#1E0200","#5B4E51","#C895C5","#320033","#FF6832","#66E1D3","#CFCDAC","#D0AC94","#7ED379","#012C58","#7A7BFF","#D68E01","#353339","#78AFA1","#FEB2C6","#75797C","#837393","#943A4D","#B5F4FF","#D2DCD5","#9556BD","#6A714A","#001325","#02525F","#0AA3F7","#E98176","#DBD5DD","#5EBCD1","#3D4F44","#7E6405","#02684E","#962B75","#8D8546","#9695C5","#E773CE","#D86A78","#3E89BE","#CA834E","#518A87","#5B113C","#55813B","#E704C4","#00005F","#A97399","#4B8160","#59738A","#FF5DA7","#F7C9BF","#643127","#513A01","#6B94AA","#51A058","#A45B02","#1D1702","#E20027","#E7AB63","#4C6001","#9C6966","#64547B","#97979E","#006A66","#391406","#F4D749","#0045D2","#006C31","#DDB6D0","#7C6571","#9FB2A4","#00D891","#15A08A","#BC65E9","#FFFFFE","#C6DC99","#203B3C","#671190","#6B3A64","#F5E1FF","#FFA0F2","#CCAA35","#374527","#8BB400","#797868","#C6005A","#3B000A","#C86240","#29607C","#402334","#7D5A44","#CCB87C","#B88183","#AA5199","#B5D6C3","#A38469","#9F94F0","#A74571","#B894A6","#71BB8C","#00B433","#789EC9","#6D80BA","#953F00","#5EFF03","#E4FFFC","#1BE177","#BCB1E5","#76912F","#003109","#0060CD","#D20096","#895563","#29201D","#5B3213","#A76F42","#89412E","#1A3A2A","#494B5A","#A88C85","#F4ABAA","#A3F3AB","#00C6C8","#EA8B66","#958A9F","#BDC9D2","#9FA064","#BE4700","#658188","#83A485","#453C23","#47675D","#3A3F00","#061203","#DFFB71","#868E7E","#98D058","#6C8F7D","#D7BFC2","#3C3E6E","#D83D66","#2F5D9B","#6C5E46","#D25B88","#5B656C","#00B57F","#545C46","#866097","#365D25","#252F99","#00CCFF","#674E60","#FC009C","#92896B")
	
	#rotation labels
	rot=90
	
	tp_x, tp_y = (0,0)
	fp_x, fp_y = (0,1)
	precsens_x, precsens_y = (0,2)
	
	prec_x, prec_y = (1,0)
	sens_x, sens_y = (1,1)
	f1score_x, f1score_y = (1,2)
	
	l1norm_x, l1norm_y = (2,0)
	avgre_x, avgre_y = (2,1)
	rrmse_x, rrmse_y = (2,2)
	
	selected_ranks = ranks.ranks
	
	all_tp = defaultdict(list)
	all_fp = defaultdict(list)
	all_precsens = defaultdict(list)
	all_sens = defaultdict(list)
	all_prec = defaultdict(list)
	all_fsco = defaultdict(list)
	all_lnor = defaultdict(list)
	all_avgre = defaultdict(list)
	all_rrmse = defaultdict(list)
	
	all = defaultdict(lambda : defaultdict(lambda: -1))
	
	for ri,r in enumerate(selected_ranks):
		rankid = ranks.getRankID(r)
		gt_profile = Profile([],[],[G.profilerank[rankid]])
		
		# Max number of negatives for the ROC curve and AUC calculation
		max_n = 0
		for tool in T:
			#If tool has the rank in the results
			profile_rank = Profile([],[],[tool.profilerank[rankid]])	
			if profile_rank.getSize():
				mn = sum([1 for name in profile_rank.getCol('TaxID') if name not in gt_profile.getCol('TaxID')])
				if mn>max_n: max_n = mn

		tool_short = []	
		for id,tool in enumerate(T):
			#If tool has the rank in the results
			profile_rank = Profile([],[],[tool.profilerank[rankid]])	
			if profile_rank.getSize():
				m = {'p':0,'n':0,'fp':[],'tp':[],'fn':0,'tpr':[],'fpr':[],'auc':0,'l1norm':0,'rrmse':0,'avgre':0}
				profile_rank.sort([('Abundance',-1)])
				m['p'] = gt_profile.getSize()
				m['n'] = max_n
				m['tp'] = np.cumsum([1 if name in gt_profile.getCol('TaxID') else 0 for name in profile_rank.getCol('TaxID')])
				m['fp'] = np.cumsum([1 if name not in gt_profile.getCol('TaxID') else 0 for name in profile_rank.getCol('TaxID')])
				m['fn'] = m['p'] - m['tp'][-1]
				m['tpr'] = [x/float(m['p']) for x in m['tp']] 
				m['fpr'] = [x/float(m['n']) for x in m['fp']]
				# Add first element to calculate AUC correctly and show first point on graph
				m['tpr'].insert(0,0)
				m['fpr'].insert(0,0)
				# Add last point to trace a line and calculate AUC correctly (Redundant when not necessary)
				m['tpr'].append(m['tpr'][-1]) # repeat last value in y axis
				m['fpr'].append(1) # 1 on x axis					
				m['auc'] = np.trapz(m['tpr'],m['fpr'])
				for gt in gt_profile.profilerank[rankid]:
					pr_ab = profile_rank.getSubSet('TaxID',gt['TaxID']).getCol('Abundance')
					pr_ab = np.sum(pr_ab) if pr_ab else 0
					m['l1norm'] += np.absolute(gt['Abundance'] - pr_ab)
					m['avgre'] += np.absolute(pr_ab - gt['Abundance'])/float(gt['Abundance'])
					m['rrmse'] += (np.absolute(pr_ab - gt['Abundance'])/float(gt['Abundance']))**2
				m['avgre'] = (1/float(m['p'])) * m['avgre']
				m['rrmse'] = math.sqrt((1/float(m['p'])) * m['rrmse'])
				
				tool_short.append(tool.ident)
				
				# rank for specific results
				if r==args.rank:
					all_tp[id].append(m['tp'][-1])
					all_fp[id].append(m['fp'][-1])
					all_precsens[id].append([m['tp'][-1]/float(m['tp'][-1]+m['fn']), m['tp'][-1]/float(m['tp'][-1]+m['fp'][-1])])

				all_prec[id].append(m['tp'][-1]/float(m['tp'][-1]+m['fp'][-1]))
				all_sens[id].append(m['tp'][-1]/float(m['tp'][-1]+m['fn']))
				all_fsco[id].append((2*m['tp'][-1])/float((2*m['tp'][-1]) + m['fp'][-1] + m['fn']))
				all_lnor[id].append(m['l1norm'])
				all_avgre[id].append(m['avgre'])
				all_rrmse[id].append(m['rrmse'])
				
				all['tp'][(id,ri)] = m['tp'][-1]
				all['fp'][(id,ri)] = m['fp'][-1]
				all['prec'][(id,ri)] = m['tp'][-1]/float(m['tp'][-1]+m['fp'][-1])
				all['sens'][(id,ri)] = m['tp'][-1]/float(m['tp'][-1]+m['fn'])
				all['fsco'][(id,ri)] = (2*m['tp'][-1])/float((2*m['tp'][-1]) + m['fp'][-1] + m['fn'])
				all['lnor'][(id,ri)] = m['l1norm']
			else:
				all_prec[id].append(None)
				all_sens[id].append(None)
				all_fsco[id].append(None)
				all_lnor[id].append(None)
				all_avgre[id].append(None)
				all_rrmse[id].append(None)

	print()
	out_arr = {}
	for pl,val in list(all.items()):
		print(pl + '\t' + '\t'.join(ranks.ranks))
		out_arr[pl] = np.zeros((len(T),len(ranks.ranks)))
		for i,tool in enumerate(T):
			sys.stdout.write(tool.ident)
			for j in range(len(ranks.ranks)):
				sys.stdout.write("\t" + str(val[i,j]))
				out_arr[pl][i,j] = val[i,j]
			print()
		print()
	np.savez(args.output_plot, tools=[t.ident for t in T], ranks=ranks.ranks, tp=out_arr['tp'], fp=out_arr['fp'], prec=out_arr['prec'], sens=out_arr['sens'], fsco=out_arr['fsco'], lnor=out_arr['lnor'])

	print("TP")
	for x,(id,val) in enumerate(all_tp.items()):
		ax[tp_x,tp_y].plot(x, val, color=colors[id], label=tool_short[id], marker=markers[id], markersize=marker_size, linestyle = 'None')
		print(tool_short[id], val)
	print()
	print("FP")
	for x,(id,val) in enumerate(all_fp.items()):
		ax[fp_x,fp_y].plot(x, val, color=colors[id], label=tool_short[id], marker=markers[id], markersize=marker_size, linestyle = 'None')
		print(tool_short[id], val)
	print()
	for id,val in list(all_precsens.items()):
		ax[precsens_x,precsens_y].plot(val[0][0], val[0][1], color=colors[id], label=tool_short[id], marker=markers[id], markersize=marker_size, linestyle = 'None')

	print("Precision")
	for id,val in list(all_prec.items()):
		ax[prec_x,prec_y].plot(list(range(len(val))), val, color=colors[id], label=tool_short[id], marker=markers[id], markersize=marker_size)
		print(tool_short[id], list(zip(ranks.ranks,val)))
	print()
	print("Sensitivity")
	for id,val in list(all_sens.items()):
		ax[sens_x,sens_y].plot(list(range(len(val))), val, color=colors[id], label=tool_short[id], marker=markers[id], markersize=marker_size)
		print(tool_short[id], list(zip(ranks.ranks, val)))
	print()
	print("F1 Score")
	for id,val in list(all_fsco.items()):
		ax[f1score_x,f1score_y].plot(list(range(len(val))), val, color=colors[id], label=tool_short[id], marker=markers[id], markersize=marker_size)
		print(tool_short[id], list(zip(ranks.ranks, val)))

	print()
	print("L1 Norm")
	for id,val in list(all_lnor.items()):
		ax[l1norm_x,l1norm_y].plot(list(range(len(val))), val, color=colors[id], label=tool_short[id], marker=markers[id], markersize=marker_size)
		print(tool_short[id], list(zip(ranks.ranks, val)))
	print()
	print("AVGRE")
	for id,val in list(all_avgre.items()):
		ax[avgre_x,avgre_y].plot(list(range(len(val))), val, color=colors[id], label=tool_short[id], marker=markers[id], markersize=marker_size)
		print(tool_short[id], list(zip(ranks.ranks, val)))
	print()
	print("RRMSE")
	for id,val in list(all_rrmse.items()):
		ax[rrmse_x,rrmse_y].plot(list(range(len(val))), val, color=colors[id], label=tool_short[id], marker=markers[id], markersize=marker_size)
		print(tool_short[id], list(zip(ranks.ranks,val)))
					
	
	margin = 0.05
	label_fontsize=12
	axis_fontsize=8
	
	#Legend
	ax[tp_x,tp_y].legend(ncol=len(T), shadow=True, fancybox=True, numpoints=1, bbox_to_anchor=(3.2, 1.4))
	
	# TP
	ax[tp_x,tp_y].set_xlabel(args.rank, rotation=rot)
	ax[tp_x,tp_y].set_xticks([])
	ax[tp_x,tp_y].set_ylabel('True Positives',fontsize=label_fontsize)
	ax[tp_x,tp_y].margins(margin)
	# FP
	ax[fp_x,fp_y].set_xlabel(args.rank, rotation=rot)
	ax[fp_x,fp_y].set_xticks([])
	ax[fp_x,fp_y].set_ylabel('False Positives',fontsize=label_fontsize)
	ax[fp_x,fp_y].margins(margin)
	# PRECISON X SENSITIVITY
	ax[precsens_x,precsens_y].plot([0,1], [0,1], ls="--")
	ax[precsens_x,precsens_y].set_xlim([0-margin,1+margin])
	ax[precsens_x,precsens_y].set_ylim([0-margin,1+margin])
	ax[precsens_x,precsens_y].set_ylabel('Precision',fontsize=label_fontsize)
	ax[precsens_x,precsens_y].set_xlabel('Sensitivity',fontsize=label_fontsize)

	
	# PRECISION
	ax[prec_x,prec_y].set_xlim([-1,len(selected_ranks)])
	ax[prec_x,prec_y].set_ylim([0-margin,1+margin])
	ax[prec_x,prec_y].set_xticks([r for r in range(len(selected_ranks))])
	ax[prec_x,prec_y].set_xticklabels(selected_ranks, rotation=rot)
	ax[prec_x,prec_y].set_ylabel('Precision',fontsize=label_fontsize)
	# SENSITIVITY
	ax[sens_x,sens_y].set_xlim([-1,len(selected_ranks)])
	ax[sens_x,sens_y].set_ylim([0-margin,1+margin])
	ax[sens_x,sens_y].set_xticks([r for r in range(len(selected_ranks))])
	ax[sens_x,sens_y].set_xticklabels(selected_ranks, rotation=rot)
	ax[sens_x,sens_y].set_ylabel('Sensitivity',fontsize=label_fontsize)	
	# F1SCORE
	ax[f1score_x,f1score_y].set_xlim([-1,len(selected_ranks)])
	ax[f1score_x,f1score_y].set_ylim([0-margin,1+margin])
	ax[f1score_x,f1score_y].set_xticks([r for r in range(len(selected_ranks))])
	ax[f1score_x,f1score_y].set_xticklabels(selected_ranks, rotation=rot)
	ax[f1score_x,f1score_y].set_ylabel('F1 Score',fontsize=label_fontsize)

	
	# L1NORM
	ax[l1norm_x,l1norm_y].set_xlim([-1,len(selected_ranks)])
	ax[l1norm_x,l1norm_y].set_xticks([r for r in range(len(selected_ranks))])
	ax[l1norm_x,l1norm_y].set_xticklabels(selected_ranks, rotation=rot)
	ax[l1norm_x,l1norm_y].set_ylabel('L1 norm',fontsize=label_fontsize)
	ax[l1norm_x,l1norm_y].margins(margin)
	# AVGRE
	ax[avgre_x,avgre_y].set_xlim([-1,len(selected_ranks)])
	ax[avgre_x,avgre_y].set_xticks([r for r in range(len(selected_ranks))])
	ax[avgre_x,avgre_y].set_xticklabels(selected_ranks, rotation=rot)
	ax[avgre_x,avgre_y].set_ylabel('AVGRE',fontsize=label_fontsize)
	ax[avgre_x,avgre_y].margins(margin)
	# RRMSE
	ax[rrmse_x,rrmse_y].set_xlim([-1,len(selected_ranks)])
	ax[rrmse_x,rrmse_y].set_xticks([r for r in range(len(selected_ranks))])
	ax[rrmse_x,rrmse_y].set_xticklabels(selected_ranks, rotation=rot)
	ax[rrmse_x,rrmse_y].set_ylabel('RRMSE',fontsize=label_fontsize)
	ax[rrmse_x,rrmse_y].margins(margin)


	fig.subplots_adjust(hspace=.3,wspace=.3)
	mpl.rcParams.update({'font.size': axis_fontsize})
	fig.set_size_inches(15,8)
	plt.savefig(args.output_plot,dpi=400)
	#plt.show()
	
if __name__ == "__main__":
	main()
