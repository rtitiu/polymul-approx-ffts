import numpy as np
from math import log2, ceil
import scipy.fft

def poly_base_decomposition(f, B): #B should be a power of two
	f = np.array(f, dtype = 'object')
	sngf = np.array(np.sign(f), dtype = 'int')
	f = abs(f)
	logB = int(log2(B))
	msk = B - 1 
	N = len(f)
	k = ceil(log2(max(f)) / log2(B))
	
	decomposed_f = np.zeros((N, k), dtype = 'int')
	for j in range(k):
		r = np.bitwise_and(f, msk)
		f = f >> logB
		r = np.where(sngf > 0, r, -r) 
		decomposed_f[:,j] = r
	return decomposed_f

def fast_mult(decomposed_f, decomposed_g, B = 2 ** 19):
	N = decomposed_f.shape[0]
	k_f = decomposed_f.shape[1]
	k_g = decomposed_g.shape[1]
	logB = int(log2(B))
	decomposed_f_fft = np.zeros((N + 1, k_f), dtype = 'complex128')
	decomposed_g_fft = np.zeros((N + 1, k_g), dtype = 'complex128')

	for i in range(k_f):
		decomposed_f_fft[:,i] = scipy.fft.rfft(decomposed_f[:, i], 2 * N)
	for i in range(k_g):
		decomposed_g_fft[:,i] = scipy.fft.rfft(decomposed_g[:, i], 2 * N)		
	
	decomposed_fg = np.zeros((2 * N, k_f + k_g - 1), dtype = 'int')	
	
	for i in range(k_f):
		for j in range(k_g):
			decomposed_fg[:,i + j] += np.round(scipy.fft.irfft(decomposed_f_fft[:,i] * decomposed_g_fft[:,j], 2 * N).real).astype(int) 
	fg_recovered = np.array([0] * (2 * N), dtype = 'object')
	for j in reversed(range(k_f + k_g - 1)):
		fg_recovered <<= logB
		fg_recovered += decomposed_fg[:,j]
	return fg_recovered	

def bd_polymul(f, g, base = 2 ** 19):
	f_dec = poly_base_decomposition(f, base)	
	g_dec = poly_base_decomposition(g, base)
	fpm = fast_mult(f_dec, g_dec, base)
	return fpm

