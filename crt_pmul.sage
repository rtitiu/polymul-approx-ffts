import numpy as np
from math import ceil, log2
import scipy.fft
from sage.arith.multi_modular import MultiModularBasis_base

def gen_prime_list(log_infinity_norm_f, log_infinity_norm_g, len_polys, B): #we need p1,p2,..,pk s.t. p1 * ..* pk > 2 * N * |f| * |g| i.e. k * log(p) ~ 2 * log_coefs + logN
	#B is a bound for the coefficients s.t. we are guaranteed that fft-complex with rounding works for N and B
	prime_list = []
	p = previous_prime(B)
	logps = log2(p)
	prime_list.append(p)
	while (logps < 1 + log2(len_polys) + log_infinity_norm_f + log_infinity_norm_g):
		p = previous_prime(p)
		prime_list.append(p)
		logps = logps + log2(p)
	return prime_list

def crt_representation(f, prime_list):
	N = len(f)
	f = np.array(f, dtype = object)
	k = len(prime_list)
	crt_f = np.zeros((N, k), dtype = object) 
	for j in range(k):
		crt_f[:,j] = f % prime_list[j]
	return crt_f	

def poly_mul_fft(crt_f, crt_g, p_list):
	N = crt_f.shape[0] #we assume both f and g have the same len = N 
	k = crt_f.shape[1]
	crt_fg = np.zeros((2 * N, k), dtype = 'complex128')

	tmp = np.zeros(N + 1, dtype = 'complex128')
	for j in range(k):
		tmp = scipy.fft.rfft(crt_f[:, j], 2 * N) * scipy.fft.rfft(crt_g[:, j], 2 * N)
		crt_fg[:,j] = scipy.fft.irfft(tmp, 2 * N)
		crt_fg[:,j] = np.round(crt_fg[:,j].real).astype(int)  % p_list[j]
	crt_fg = np.array(crt_fg, dtype = int)	
	return crt_fg

def crt_with_precomputation(mmb, poly_in_crt): 
	f = []
	for i in range(len(poly_in_crt)):
		f.append(mmb.crt(poly_in_crt[i,:])) 
	return np.array(f) 	

def crt_polymul(f,g,prime_list, mmb):
	crt_f = crt_representation(f, prime_list)
	crt_g = crt_representation(g, prime_list)
	crt_fg = poly_mul_fft(crt_f, crt_g, p_list)
	fg = crt_with_precomputation(mmb, crt_fg)
	return np.array(fg)	