from dolfin import *
import numpy
import nose.tools as nt
from scitools.std import array

def main(nx, dt):
	# Parameters
	T = 1.0
	rho = 1.0
	# Create mesh and define function space
	mesh = UnitSquareMesh(nx, nx)
	V = FunctionSpace(mesh, 'Lagrange', 1)
	
	# Define boundary conditions
	u0 = Expression('exp(-pi*pi*t)*cos(pi*x[0])', t=0)
	# Initial condition
	u_1 = interpolate(u0, V)
	
	
	#Define variational problem
	u = TrialFunction(V)
	v = TestFunction(V)
	alpha = Expression('1')
	f = 0
	a = u*v*dx + (dt/rho)*alpha*inner(nabla_grad(u), nabla_grad(v))*dx
	L = u_1*v*dx + (dt/rho)*f*v*dx
	
	A = assemble(a)
	b = None
	
	u = Function(V)
	
	# Initial condition already stored in u_1, so start iteration at t=dt
	t = dt
	E = [];
	while t<=T:
		b = assemble(L, tensor=b)
		solve(A, u.vector(), b)
		u0.t = t
		u_e = interpolate(u0, V)	
		e = u_e.vector().array() - u.vector().array()
		E.append(numpy.sqrt(numpy.sum(e**2)/u.vector().array().size))
		#plot(u_1)
		u_1.assign(u) # Copy solution to u_1 to prepare for next time-step
		t += dt
	return array(E).max()
	
	
def test_conv():
	E = []
	h = []
	r = []
	for nx in [4, 8, 16, 32, 64]:
		print ('Computing for nx=%s' % nx)
		h.append(1.0/nx**2)
		dt = h[-1]
		E.append(main(nx, dt))
		
	print numpy.array(E)/numpy.array(h)
	nt.assert_almost_equal(E[-1]-E[-2], 0, places=2)
	for i in range(len(E)-1):
		R = ln(E[i]/E[i-1])/ln(h[i]/h[i-1])
		r.append(R)
	print r
	
	
	
if __name__ == '__main__':
	set_log_level(50)
	test_conv()
