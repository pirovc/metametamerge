from metametamerge.Profile import Profile
import numpy as np

class Tools(Profile):

	def __init__(self, file, ident, method, profile, ranks, db=None):	
		self.file = file
		self.ident = ident
		self.method = method
		Profile.__init__(self, profile, ranks)
		if db: 
			#check presence on it's on database profile
			self.checkDB(db, ranks)
			# Estimate abundance for binning methods
			if self.method=='b': self.estimateAbundance(db)
		# Normalize abundance for all
		self.normalizeAbundance()

	def checkDB(self, db, ranks):
		db_taxids = db.getCol('TaxID')
		taxids_not_found = []
		for profilerank in list(self.profilerank.values()): 
			taxids = profilerank.getCol('TaxID')
			taxids_in_db = np.in1d(taxids,db_taxids)
			for t in taxids[~taxids_in_db]:
				print(("tools TaxID not found in the database - ignoring entry [%d] rank [%s]" % (t,ranks.getRankName(profilerank.rankid))))
			# Filter them out
			profilerank.filter(taxids_in_db)
			
	def estimateAbundance(self,db): 
		for profilerank in list(self.profilerank.values()): 
			for pr in profilerank:
				taxid = pr['TaxID']
				old_ab = pr['Abundance']
				db_taxid = db.getSubSet('TaxID',pr['TaxID'])

				# Sum duplicated taxids in database cause by inconsistencies in the NCBI databases (usually old accessions) providing different tax. names so they cannot be merged
				if db_taxid.getSize()>1: print("tools estimateAbundance WARNING taxid duplicated in the database ", taxid)
				taxid_db_len = np.sum(db_taxid.getCol('Abundance'))
				new_ab = old_ab/float(taxid_db_len)

				profilerank.profilerank[profilerank.getCol('TaxID')==pr['TaxID'],profilerank.cols['Abundance']] = new_ab
			
	def normalizeAbundance(self): 
		for profilerank in list(self.profilerank.values()): 
			totalAbundance = np.sum(profilerank.getCol('Abundance'))
			# Normalize by the max (everything always between 0 and 1)
			if totalAbundance: profilerank.profilerank[:,profilerank.cols['Abundance']] *= 1/float(totalAbundance)
			
				
	def filterMaxResults(self, max_results):
		for profilerank in list(self.profilerank.values()): 
			if profilerank.getSize():
				# Sort by abundance
				profilerank.sort([('Abundance',-1)])
				# Select max number number of results
				mr = profilerank.getSize() if profilerank.getSize()<max_results else max_results
				max_idx = np.array([True]*mr + [False]*(profilerank.getSize()-mr))
				# Filter results
				profilerank.filter(max_idx)
	
	def filterMinRelativeAbundance(self, min_relative_abundance):
		for profilerank in list(self.profilerank.values()): 
			min_ra_idx = profilerank.getCol('Abundance')>=min_relative_abundance
			# Filter results
			profilerank.filter(min_ra_idx)
	
	def filterRepeatedNameRank(self):
		from collections import Counter
		for profilerank in list(self.profilerank.values()): 
			repeated_taxids = [item for item, count in Counter(profilerank.getCol('TaxID')).items() if count > 1]
			for r in repeated_taxids:
				print("tools filterRepeatedNameRank WARNING ignorning repeated taxid in same rank", r, self.file)
				repeated_taxids_idx = profilerank.getCol('TaxID')==r
				profilerank.filter(~repeated_taxids_idx)
	
