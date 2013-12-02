from dolfin import *
import numpy
import nose.tools as nt
import time


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
	
	f = Expression('-rho*pow(x[0],3.)/3. + rho*pow(x[0],2.)/2. +\
		8*pow(t,3.)*pow(x[0],7.)/9. - 28*pow(t,3.)*pow(x[0],6.)/9. +\
		7*pow(t,3)*pow(x[0],5.)/2. -5*pow(t,3.)*pow(x[0],4.)/4. + 2*t*x[0] - t',
		rho=rho, t=0)
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
		E1 = numpy.sqrt(numpy.sum(e**2)/u.vector().array().size)
		if(counter%10 == 0):
			E.append(E1)
			#plot(u)
		u_1.assign(u) # Copy solution to u_1 to prepare for next time-step
		t += dt
		counter += 1
	return E
	
def test_compare():
	E = main(20, 0.001)
	for i in range(len(E)):
		nt.assert_almost_equal(E[i], 0, places = 2)
	
		
		
if __name__ == '__main__':
	set_log_level(50)
	test_compare()
