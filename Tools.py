#from metametamerge.Profile import Profile
from Profile import Profile
import numpy as np

class Tools(Profile):
		
	def __init__(self, file, ident, method, profile, ranks, verbose, db=None):	
		self.file = file
		self.ident = ident
		self.method = method
		Profile.__init__(self, profile, ranks)
		if db: 
			#check presence on it's on database profile
			self.checkDB(db, ranks, verbose)
			# Estimate abundance for binning methods
			if self.method=='b': self.estimateAbundance(db, verbose)
		# Normalize abundance for all
		self.normalizeAbundance()

	def checkDB(self, db, ranks, verbose):
		db_taxids = db.getCol('TaxID')
		taxids_not_found = []
		for profilerank in list(self.profilerank.values()): 
			taxids = profilerank.getCol('TaxID')
			taxids_in_db = np.in1d(taxids,db_taxids)
			if verbose:
				for t in taxids[~taxids_in_db]: print(("Ignored entry [%d] rank [%s] - taxon not found in the database profile" % (t,ranks.getRankName(profilerank.rankid))))
			if taxids[~taxids_in_db].any(): print(("\t%d filtered taxons [%s] (not found on database profile)") % (len(taxids[~taxids_in_db]),ranks.getRankName(profilerank.rankid)))
			# Filter them out
			profilerank.filter(taxids_in_db)
			
	def estimateAbundance(self, db, verbose): 
		for profilerank in list(self.profilerank.values()): 
			for pr in profilerank:
				taxid = pr['TaxID']
				# Get db length for taxid (unique after mergeRepeatedTaxIDs)
				db_len = db.getSubSet('TaxID',taxid).getCol('Abundance')[0]
				new_ab = pr['Abundance']/float(db_len)
				profilerank.profilerank[profilerank.getCol('TaxID')==taxid,profilerank.cols['Abundance']] = new_ab
			
	def normalizeAbundance(self): 
		for profilerank in list(self.profilerank.values()): 
			totalAbundance = np.sum(profilerank.getCol('Abundance'))
			# Normalize by the max (everything always between 0 and 1)
			if totalAbundance: profilerank.profilerank[:,profilerank.cols['Abundance']] *= 1/float(totalAbundance)
			
				
	def filterMaxResults(self, max_results, ranks):
		for profilerank in list(self.profilerank.values()): 
			if profilerank.getSize():
				# Sort by abundance
				profilerank.sort([('Abundance',-1)])
				# Select max number number of results
				mr = profilerank.getSize() if profilerank.getSize()<max_results else max_results
				max_idx = np.array([True]*mr + [False]*(profilerank.getSize()-mr))
				print(("\t%s - %d entries filtered [%s]") % (self.ident, profilerank.profilerank[~max_idx].shape[0], ranks.getRankName(profilerank.rankid)))
				# Filter results
				profilerank.filter(max_idx)
	
	def filterMinRelativeAbundance(self, min_relative_abundance, ranks):
		for profilerank in list(self.profilerank.values()): 
			if profilerank.getSize():
				min_ra_idx = profilerank.getCol('Abundance')>=min_relative_abundance
				print(("\t%s - %d entries filtered [%s]") % (self.ident, profilerank.profilerank[~min_ra_idx].shape[0], ranks.getRankName(profilerank.rankid)))
				# Filter results
				profilerank.filter(min_ra_idx)
