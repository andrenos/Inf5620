import sys
from scitools.std import zeros, linspace, sin, cos ,tan, sqrt, log
import matplotlib.pyplot as plt

def define_command_line_options():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--T', '--stop_time', type=float, 
                    default=20.0, help='end time of simulation',
						metavar='t')
	parser.add_argument('--dt', type=float,
    					default=0.1, help='timestep for the discrete apporoximation',
    					metavar='dt')
	parser.add_argument('--v0', '--initial_condition', type=float,
    					default=-0.0, help='initial condition v(0)',
    					metavar='v0')
	parser.add_argument('--makeplot', action='store_true',
    					help='display plot or not')
	parser.add_argument('--a', type=float,
    					default=0.003, help='coefficient in ODE',
    					metavar='a')
	parser.add_argument('--m', '--body_mass', type=float,
    					default=100., help='body mass', metavar='m')
	return parser


def skydiving_iterate(v, t, dt, X, Y):
    """
    Problem specific function. Implements the needed computation for 
    doing one time-step iteration.

    """
    return (v + dt*X(t))/(1 + dt*Y(t)*abs(v))

def solver(T, dt, v0, a, m, Source=None):
    v = zeros(int(T/dt) + 1)
    v[0] = v0; g = 9.81;
    def X(t):
        if(Source == None):
            return -g
        else:
            return -g + Source(t+dt/2)/m
            #return -g + 0.5*(Source(t) + Source(t+dt))/m
        
    def Y(t):
        return a
    
    for i in range(1,len(v)):
        v[i] = skydiving_iterate(v[i-1], dt*(i-1), dt, X, Y)
    return v, linspace(0, T, T/dt +1)


def plot_v(t, v):
    p1 = plt.plot(t,v)
    plt.xlabel('Time [s]')
    plt.ylabel('Velocity [m/s]')
    plt.title('Velocity for the skydiver as a function of time')
    plt.show()

def plot_forces(t, v, m, a):
    plt.figure()
    drag = -m*a*abs(v)*v
    grav = [-m*9.81]*len(v)
    Boyancy = [1. * 9.81 * 0.1]*len(v) # rho * g * V
    Fsum = drag+grav+Boyancy
    plt.plot(t, drag, t, grav, t, Boyancy,  t, Fsum)
    plt.legend(["Drag force", "Gravity force", "Boyancy", "Sum of forces"])

def main():
    parser = define_command_line_options()
    args = parser.parse_args()
        
    v, t = solver(args.T, args.dt, args.v0, args.a, args.m)
    for i in range(len(v)):
        print t[i], v[i]
    if(args.makeplot):
        plot_v(t, v)
        plot_forces(t, v, args.m, args.a)
        raw_input()

if(__name__ == "__main__"):
	main()

def test_constant():
    """
    Test for the solvers ability to reproduce a constant solution
    """
    import nose.tools as nt
    A = -0.11; g = 9.81; m = 50.; a = 1.0; T = 10.; dt = 0.01;
    def src(t):
        return m*g + m*a*abs(A)*A
    v, t = solver(T, dt, A, a, m, Source=src)
    for i in range(len(v)):
        nt.assert_almost_equal(A, v[i], delta=1e-12)
        
def test_linear():
    """
    Test for the solvers ability to reproduce a linear solution
    by fitting the source term so that a linear solution should indeed
    be a solution to the differential equation
    """
    import nose.tools as nt
    A = -0.11; B = -0.13; g = 9.81; m = 50.; a = 0.91; T = 10.; dt = 0.01;

    def exact(t):
        return A*t+B

    def src(t):
        return m*g + m*a*abs(exact(t-dt/2.))*exact(t+dt/2.) + m*A
        
    v, t = solver(T, dt, B, a, m, Source=src)
    ve = exact(t)
    diff = abs(ve - v)
    nt.assert_almost_equal(diff.max(), 0, delta=1e-12)
        
def test_terminalVelocity():
    """
    Test for the terminal velocity with no source term.
    """
    import nose.tools as nt
    T = 10.; dt = 0.1; g = 9.81; m = 50.; a = 1.0;
    v, t = solver(T, dt, -0.1, a, m)
    nt.assert_almost_equal(v[-1], -sqrt(g/a), delta=1e-8)
        
def test_convergenceRates():
    dt_start = 1.0; num_dt = 10
    E_values = zeros(num_dt)
    T = 10.; g = 9.81; m = 50.; a = 0.91;
    dt = zeros(num_dt); dt[0] = dt_start
    for i in range(1,len(dt)):
        dt[i] = dt[i-1]/2.
    print "dt=", dt
    
    A = -0.39; B = 0.76; C = -0.145
    def exact(t):
        return A*t**3  + B*t + C
    def src(t):
        return m*g + m*a*abs(exact(t))*exact(t) + m*(3*A*t**2 + B)
    
    for i in range(num_dt):
        v, t = solver(T, dt[i], exact(0), a, m, src)
        ve = exact(t)
        diff = v - ve
        E_values[i] = sqrt(dt[i]*sum(diff**2))

    r=zeros(len(E_values)-1)
    for i in range(1, len(r)):
        r[i] = (log(E_values[i-1]/E_values[i]))/(log(dt[i-1]/dt[i]))
    print("E=", E_values)
    print("r=", r)
    import nose.tools as nt
    nt.assert_almost_equal(r[-1], 2, delta=0.1)