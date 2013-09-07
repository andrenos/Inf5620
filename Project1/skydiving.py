import sys
from scitools.std import zeros, linspace, sqrt, log, concatenate
import matplotlib.pyplot as plt

def define_command_line_options():
    """ 
    Function for giving command line argument as option-value pairs
    Usage: >>python skydiving.py --T 60 --v0 0.
    etc.
    """
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--T', '--stop_time', type=float, 
                        default=20.0, help='end time of simulation', 
                        metavar='t')
    parser.add_argument('--dt', type=float, default=0.1,
                        help='timestep for the discrete apporoximation',
                        metavar='dt')
    parser.add_argument('--v0', '--initial_condition', type=float,
                        default=-0.0, help='initial condition v(0)',
                        metavar='v0')
    parser.add_argument('--makeplot', action='store_true',
                        help='display plot or not')
    parser.add_argument('--rho', type=float, default=1.0,
                        help='air mass density', metavar='rho')
    parser.add_argument('--Cd', type=float, default=1.2,
                        help='drag coefficient', metavar='Cd')
    parser.add_argument('--m', '--body_mass', type=float, default=100.,
                        help='body mass', metavar='m')
    parser.add_argument('--A', type=float, default=0.5,
                        help='body cross sectional area',
                        metavar='A')
    parser.add_argument('--tp', type=float, default=-1,
                        help='time of parachute deployment', metavar='tp')
    return parser
    
    
def plot_v(t, v):
    """
    Function for plotting the velocity of the skydiver as a function of time
    """
    p1 = plt.plot(t,v)
    plt.xlabel('Time [s]')
    plt.ylabel('Velocity [m/s]')
    plt.title('Velocity for the skydiver as a function of time')
    plt.show()
    plt.savefig('Parachute_velocity.png')

def plot_forces(t, v, m, a):
    """
    Function for plotting the forces when running 
    the program without parachute deployment
    """
    plt.figure()
    drag = -m*a*abs(v)*v
    grav = [-m*9.81]*len(v)
    Boyancy = [1. * 9.81 * 0.1]*len(v) # rho * g * V
    Fsum = drag+grav+Boyancy
    plt.plot(t, drag, t, grav, t, Boyancy,  t, Fsum)
    plt.legend(["Drag force", "Gravity force", "Boyancy", "Sum of forces"])
    plt.savefig('Forces.png')

def plot_forces_parachute(t, v, dt, tp, m, a_first, a_last):
    """
    Function for plotting the forces when running 
    the program with parachute deployment
    """
    plt.figure()
    drag = zeros(len(v))
    for i in range(len(v)):
        if(i*dt <= tp):
            drag[i] = -m*a_first*abs(v[i])*v[i]
        else:
            drag[i] = -m*a_last*abs(v[i])*v[i]
    grav = [-m*9.81]*len(v)
    Boyancy = [1. * 9.81 * 0.1]*len(v) # rho * g * V
    Fsum = drag+grav+Boyancy
    plt.plot(t, drag, t, grav, t, Boyancy,  t, Fsum)
    plt.legend(["Drag force", "Gravity force", "Boyancy", "Sum of forces"])
    plt.savefig('Parachute_forces.png')

def plot_x(t, x):
    """
    Function for plotting the vertical position of the skydiver
    as a function of time
    """
    plt.figure()
    plt.plot(t, x)
    plt.title("Vertical position of the skydiver as a function of time")
    plt.xlabel("Time t [s]")
    plt.ylabel("Height [m]")
    plt.savefig('Parachute_position.png')


def skydiving_iterate(v, t, dt, X, Y):
    """
    Problem specific function. Implements the needed computation for 
    doing one time-step iteration. Called from solver and solver_parachute.

    """
    return (v + dt*X(t))/(1 + dt*Y(t)*abs(v))

def solver(T, dt, v0, Cd, rho, A, m, Source=None):
    """
    This is the main solver function for the program. It takes the problem 
    specific variables as arguments, and returns two meshes containing the 
    velocity and time points respectively.
    """
    a = Cd*rho*A/(2*m) # Set up the constant a for compact code
    v = zeros(int(T/dt) + 1) # Create the velocity mesh
    v[0] = v0; g = 9.81;#Initial velocity and the value of gravity acceleration
    
    # Description of the functions X(t) and Y(t) is given in the PDF file
    def X(t):  
        if(Source == None):
            return -g
        else:
            return -g + Source(t+dt/2.)/m
        
    def Y(t):
        return a
    
    #Calculate the velocity at each meshpoint
    for i in range(1,len(v)):
        v[i] = skydiving_iterate(v[i-1], dt*(i-1), dt, X, Y)
    return v, linspace(0, T, T/dt +1)

def solver_parachute(T, dt, v0, tp, Cd, rho, A, m):
    """
    This is an extension of the main solver function that includes the 
    deployment of the parachute. It requires a separate implementation 
    because the final time is not known, and the constants Cd and A are
    time-dependent. It also returns the position mesh
    """
    
    # Mesh initialization
    v = []; v.append(v0); t=[]; t.append(0);
    x = []; x.append(3000.)
    
    def X(t):
        return -9.81
    def Y(t):
        return Cd*rho*A/(2.*m)
        
    # Simulate the same way as solver() untill time tp
    while t[-1]<=tp:
        v.append(skydiving_iterate(v[-1], t[-1], dt, X, Y))
        x.append(x[-1] + dt*v[-1])
        t.append(t[-1] + dt)
        
    # Z is the same as the function Y, but with different values
    def Z(t):
        return 1.8*1.0*44./(2.*m)
    
    # Continue the simulation untill the skydiver hits the ground
    while x[-1] > 0:
        v.append(skydiving_iterate(v[-1], t[-1], dt, X, Z))
        x.append(x[-1] + dt*v[-1])
        t.append(t[-1] + dt)
    return v, t, x


def main():
    parser = define_command_line_options()
    args = parser.parse_args()
    
    # Check if the parachute should be deployed or not
    if(args.tp == -1):
        # No parachute deployment
        v, t = solver(args.T, args.dt, args.v0, args.Cd, args.rho, args.A,
                    args.m)
        if(args.makeplot):
            plot_v(t, v)
            plot_forces(t, v, args.m, args.Cd*args.rho*args.A/(2.*args.m))
            raw_input()
    else:
        v, t, x = solver_parachute(args.T, args.dt, args.v0, args.tp, args.Cd,
                                args.rho, args.A, args.m)
        if(args.makeplot):
            plot_v(t, v)
            plot_forces_parachute(t, v, args.dt, args.tp, args.m, 
                                args.Cd*args.rho*args.A/(2.*args.m),
                                1.8*args.rho*44/(2.*args.m))
            plot_x(t, x)
            raw_input()

if(__name__ == "__main__"):
	main()

def test_constant():
    """
    Test for the solvers ability to reproduce a constant solution.
    """
    import nose.tools as nt
    
    C = -0.11; g = 9.81; m = 50.; T = 10.; dt = 0.01;
    Cd = 1.2; rho = 1.0; A = 0.5; 
    a = Cd*rho*A/(2.*m)
    def src(t):
        return m*g + m*a*abs(C)*C
    v, t = solver(T, dt, C, Cd, rho, A, m, Source=src)
    for i in range(len(v)):
        nt.assert_almost_equal(C, v[i], delta=1e-12)
        
def test_linear():
    """
    Test for the solvers ability to reproduce a linear solution
    by fitting the source term so that a linear solution should indeed
    be a solution to the differential equation
    """
    import nose.tools as nt
    A = -0.11; B = -0.13; g = 9.81; m = 50.; T = 10.; dt = 0.01;
    Cd = 1.2; rho = 1.0; A = 0.5;
    a = Cd*rho*A/(2.*m)
    def exact(t):
        return A*t+B

    def src(t):
        return m*g + m*a*abs(exact(t-dt/2.))*exact(t+dt/2.) + m*A
        
    v, t = solver(T, dt, B, Cd, rho, A, m, Source=src)
    ve = exact(t)
    diff = abs(ve - v)
    nt.assert_almost_equal(diff.max(), 0, delta=1e-12)
        
def test_terminalVelocity():
    """
    Test for the terminal velocity with no source term.
    """
    import nose.tools as nt
    T = 30.; dt = 0.1; g = 9.81; m = 50.;
    Cd = 1.2; rho = 1.0; A = 0.5;
    a = Cd*rho*A/(2.*m)
    v, t = solver(T, dt, -0.1, Cd, rho, A, m)
    nt.assert_almost_equal(v[-1], -sqrt(g/a), delta=1e-4)
        
def test_convergenceRates():
    """
    Test for the convergence rates of the solver.
    The expected result is that the variable r takes the value 2, because
    the Crank-Nicolson scheme and the geometric average have errors of order
    dt**2. The final error should then be O(dt**2) which gives r=2. 
    """
    dt_start = 1.0; num_dt = 10
    E_values = zeros(num_dt)
    T = 10.; g = 9.81; m = 50.; Cd = 1.2; rho = 1.0; A = 0.5;
    a = Cd*rho*A/(2.*m)
    
    dt = zeros(num_dt); dt[0] = dt_start
    for i in range(1,len(dt)):
        dt[i] = dt[i-1]/2.
    
    D = -0.39; B = 0.76; C = -0.145
    def exact(t):
        return D*t**3  + B*t + C
    def src(t):
        return m*g + m*a*abs(exact(t))*exact(t) + m*(3*D*t**2 + B)
    
     # Simulate for different timesteps, and store the error in E_values
    for i in range(num_dt):
        v, t = solver(T, dt[i], exact(0), Cd, rho, A, m, src)
        ve = exact(t)
        diff = v - ve
        E_values[i] = sqrt(dt[i]*sum(diff**2))
    
    # Calculate r-values corresponding to the error with each different timestep
    r=zeros(len(E_values)-1)
    for i in range(1, len(r)):
        r[i] = (log(E_values[i-1]/E_values[i]))/(log(dt[i-1]/dt[i]))
    import nose.tools as nt
    nt.assert_almost_equal(r[-1], 2, delta=0.1)
    

    