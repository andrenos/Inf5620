from scitools.std import *

	
def cooling(T0, k, T_s, t_end, dt, theta=0.5):
	Nt = int(t_end / dt)
	T = zeros((Nt, len(T0)))
	for i in range(len(T0)):
		T[0][i] = T0[i];
		
	for i in xrange(Nt-1):
		for j in range(len(T0)):
			nom = T[i][j]*(-k*theta + 1./dt) + T_s*k
			denom = 1./dt + (1 - theta)*k
			T[i+1][j] = nom/denom
	return T
	
if(__name__ == "__main__"):
	T = cooling([200, 300], 1, 50, 10, 1, 0.5);
	for i in range(len(T)):
		print T[i][0], T[i][1]
	
			