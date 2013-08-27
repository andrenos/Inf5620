import sys
from scitools.std import *

def skydiving_iterate(v, dt, k, m, source):
	g = 9.81
	sourceNom = 0.
	sourceDenom = 0.
	nom = -g + sourceNom + v/dt
	denom = 1./dt+ k*abs(v) + sourceDenom;
	return nom/denom
	
def solver(t, dt, v0, k, m, source):
	v = zeros(int(t/dt) + 1)
	v[0] = v0
	for i in range(1,int(t/dt)+1):
		v[i] = skydiving_iterate(v[i-1], dt, k, m, source)
	return v, linspace(0, t, t/dt +1)


if(__name__ == "__main__"):
	v, t = solver(float(sys.argv[1]), float(sys.argv[2]),
				float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), 0.)
	for i in range(len(v)):
		print t[i], v[i]

