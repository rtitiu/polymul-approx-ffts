load("crt_pmul.sage")
from bd_pmul import bd_polymul
from time import time
from random import SystemRandom
from flint import fmpz_poly


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

log_coefs = 200
bound = 2 ^ (log_coefs)

p_list = gen_prime_list(log_coefs, log_coefs, N, B)

mm = MultiModularBasis_base(p_list) 
time()
time()
f = uniform_vector(-bound, bound, N)
g = uniform_vector(-bound, bound, N)

ti_bd = time()
fg_flint = flint_polymul(f,g)
tf_bd = time()

ti_crt = time()
fg_crt = crt_polymul(f,g, p_list, mm)
tf_crt = time()

nr_zros = len(fg_crt) - len(fg_flint)
fg_flint = np.append(fg_flint, [0] * nr_zros)

fg_flint = [int(_) for _ in fg_flint]
fg_crt = [int(_) for _ in fg_crt]

if not np.array_equal(fg_flint, fg_crt):
	print("Polynomial multiplication not correct!")
else:
	print("N = {}, k = {} \n".format(N, log_coefs))		
	print("flint polymul time : {:.3f}".format(tf_bd - ti_bd))
	print("crt_polymul time   : {:.3f}".format(tf_crt - ti_crt))
