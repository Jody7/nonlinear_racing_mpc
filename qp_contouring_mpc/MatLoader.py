from scipy.io import loadmat
import scipy.io as spio

class MatLoader():
	def loadmat(self, filename):
		data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
		return self._check_keys(data)

	def _check_keys(self, dict):
		for key in dict:
			if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
				dict[key] = self._todict(dict[key])
		return dict        

	def _todict(self, matobj):
		dict = {}
		for strg in matobj._fieldnames:
			elem = matobj.__dict__[strg]
			if isinstance(elem, spio.matlab.mio5_params.mat_struct):
				dict[strg] = self._todict(elem)
			else:
				dict[strg] = elem
		return dict