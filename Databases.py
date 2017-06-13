#from metametamerge.Profile import Profile
from Profile import Profile

class Databases(Profile):

	def __init__(self, file, profile, ranks):
		self.file = file
		Profile.__init__(self, profile, ranks)
		
		# Merge repeated taxids (when taxid changes in the new taxonomy version)
		merged_taxids = self.mergeRepeatedTaxIDs()
		if merged_taxids: print(("\t%d taxons with merged entries [%s]") % (len(merged_taxids),merged_taxids))
			
