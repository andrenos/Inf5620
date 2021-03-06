from dolfin import *
import numpy
import nose.tools as nt
import time
from scitools.std import array
import scitools.std as sc

def main(nx, dt):
	# Parameters
	T = 1.0
	rho = 1.0
	# Create mesh and define function space
	mesh = UnitIntervalMesh(nx)
	V = FunctionSpace(mesh, 'Lagrange', 1)
	
	# Define boundary conditions
	u0 = Expression('t * x[0]*x[0]*(0.5 - x[0]/3.)', t=0)
	# Initial condition
	u_1 = interpolate(u0, V)
	
	
	#Define variational problem
	u = TrialFunction(V)
	v = TestFunction(V)
	
	f = Expression('rho*x[0]*x[0]*(-2*x[0] + 3)/6. - \
		(-12*t*x[0] + 3*t*(-2*x[0] + 3))*(pow(x[0],4.)*\
		pow((-dt + t),2)*pow((-2*x[0] + 3),2.) + 36)/324.\
		- (-6*t*x[0]*x[0] + 6*t*x[0]*(-2*x[0] + 3))*\
		(36*pow(x[0],4.)*pow((-dt + t),2.)*(2*x[0] - 3)+\
		36*pow(x[0],3.)*pow((-dt + t),2.)*pow((-2*x[0] + 3),2.))/5832.',
		rho=rho, dt=dt, t=0)
	a = u*v*dx + (dt/rho)*(1-u_1**2)*inner(nabla_grad(u), nabla_grad(v))*dx		
	L = u_1*v*dx + (dt/rho)*f*v*dx 
	
	u = Function(V)
	
	# Initial condition already stored in u_1, so start iteration at t=dt
	t = dt
	counter = 1
	E=[]
	while t<=T:
		u0.t = t; f.t = t
		solve(a==L, u)
		u_e = interpolate(u0, V)
		e = u_e.vector().array() - u.vector().array()
		E.append(numpy.sqrt(numpy.sum(e**2)/u.vector().array().size))
		if(counter%10 == 0):
			print e
		
		u_1.assign(u) # Copy solution to u_1 to prepare for next time-step
		t += dt
		counter += 1
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
	nt.assert_almost_equal(E[-1]-E[-2], 0, places=4)
	for i in range(len(E)-1):
		R = ln(E[i]/E[i-1])/ln(h[i]/h[i-1])
		r.append(R)
	print r
	
	
	
	
if __name__ == '__main__':
	set_log_level(50)
	test_conv()
