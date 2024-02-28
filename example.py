from bd_pmul import bd_polymul
from random import SystemRandom
from time import time
from flint import fmpz_poly
import numpy as np

def uniform_vector(A,B,N): #returns a vector of len N with uniform entries in [A, B]  
	return [int(SystemRandom().randrange(B - A + 1) + A) for _ in range(N)]

def flint_polymul(a,b):
	f = fmpz_poly(a)
	g = fmpz_poly(b)
	return f * g

NN = [1024,2048, 4096, 8192]
logBB = [20,20,19,18]

N = NN[1]
logB = logBB[1]
B = 2 ** logB

k = 200
coefficient_bound = 2 ** k

time()

f = []
g = []

f = uniform_vector(- coefficient_bound, coefficient_bound, N)
g = uniform_vector(- coefficient_bound, coefficient_bound, N)

ti_flint = time()
fg_flint = flint_polymul(f,g)
tf_flint = time()

ti_bd = time()
fg = bd_polymul(f,g,B)
tf_bd = time()

fg_flint = np.array(fg_flint)
nr_zros = len(fg) - len(fg_flint)
fg_flint = np.append(fg_flint, [0] * nr_zros)

if not np.array_equal(fg_flint, fg):
	print("\n Polynomial multiplication is not correct!")
else:
	print("\nN = {}, k = {}\n".format(N,k))	
	print("flint time     : {:.3f}s".format(tf_flint - ti_flint))
	print("bd_polymul time: {:.3f}s".format(tf_bd - ti_bd))