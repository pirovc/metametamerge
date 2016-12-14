class Ranks:

	def __init__(self,ranks):
		self.ranks = ranks
		
	def getRankID(self,rank):
		return self.ranks.index(rank)
	
	def getIDs(self):
		return list(range(len(self.ranks)))
		
	def getRankName(self,rankid):
		return self.ranks[rankid]
