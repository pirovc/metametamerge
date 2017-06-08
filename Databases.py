#from metametamerge.Profile import Profile
from Profile import Profile

class Databases(Profile):

	def __init__(self, file, profile, ranks):
		self.file = file
		Profile.__init__(self, profile, ranks)
