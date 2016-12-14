import numpy as np

class Profile:
	def __init__(self,profile,ranks,profileRanks=None):
		self.profilerank = {}

		if profileRanks!=None:
			# from getSubSet (array of ProfileRank)
			for pr in profileRanks:
				self.profilerank[pr.rankid] = pr
		else:
			# From parse_files
			for rankid in ranks.getIDs():
				#Extract only iterated rankid and other columns that are not the rank column (redundant)
				self.profilerank[rankid] = ProfileRank(profile[np.ix_(profile[:,1]==rankid, [0,2,3])],rankid,sum(profile[:,1]==rankid))

	def __iter__(self):
		for rankid,profilerank in list(self.profilerank.items()):
			yield rankid,profilerank
			
	def getSize(self,rankid=None):
		if rankid!=None:
			return self.profilerank[rankid].getSize()
		else:
			return int(np.sum([profilerank.getSize() for profilerank in list(self.profilerank.values())]))
	
	def getSubSet(self,col,val):
		subset_pr = []
		for profilerank in list(self.profilerank.values()): 
			r = profilerank.getSubSet(profilerank.getCol(col)==val)
			if r.getSize():	subset_pr.append(r)
		return Profile([],[],subset_pr)

	def getCol(self,col):
		vals = []
		for profilerank in list(self.profilerank.values()): 
			vals.extend(profilerank.getCol(col))
		return vals
	
	def sort(self,cols):
		for profilerank in list(self.profilerank.values()):
			profilerank.sort(cols)

		
class ProfileRank:
	columns = ['Presence','TaxID','Abundance']
	def __init__(self,profileRank,rankid,original_len=None):
		self.cols = {c:i for i,c in enumerate(ProfileRank.columns)}
		self.rankid = rankid
		self.profilerank = profileRank
		if original_len: self.original_len = original_len

	def __iter__(self):
		for v in self.profilerank:
			yield  {c:v[i] for i,c in enumerate(self.columns)}
			
	def getSize(self):
		return int(self.profilerank.shape[0])
		
	def getSubSet(self,idx):
		return ProfileRank(self.profilerank[idx],self.rankid)
	
	def filter(self,idx):
		self.profilerank = self.profilerank[idx]
	
	def getCol(self,col):
		return self.profilerank[:,self.cols[col]]
		
	def getRow(self,start,end):
		return ProfileRank(self.profilerank[start:end,:],self.rankid)
				
	def sort(self,cols):
		ord = tuple(self.profilerank[:,self.cols[c]]*s for c,s in cols)
		sort_idx = np.lexsort(ord)
		self.profilerank = self.profilerank[sort_idx]
