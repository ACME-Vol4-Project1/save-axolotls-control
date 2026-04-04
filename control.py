# TODO: add control utilities

import numpy as np
import sympy as sy
from naming import *                            # For the model definitions
from model import sy_params_dynamic, sy_f_simple
# from model import sy_f_full as f              # Defining them both here for the time being

from scipy.integrate import solve_ivp, solve_bvp
from scipy import linalg as la
import model as m
from matplotlib import pyplot as plt

### FOR SOLVE_BVP ###
def define_variables(no_dead=True):
    """Define sympy variables for use in other functions."""
    u1, u2 = sy.symbols("u1, u2")
    S, I, T1, D = sy.symbols("S, I, T1, D")
    u = sy.Matrix([u1, u2])

    if no_dead:
        x = sy.Matrix([S, I, T1])
        lam = sy.symbols("l1, l2, l3")
    else:
        x = sy.Matrix([S, I, T1, D])
        lam = sy.symbols("l1, l2, l3, l4")
    return x, u, sy.Matrix([lam])

def L(x, u, no_dead=True, u1="const"):
    """Inside of the cost integral, x.TQx + u.TRu."""
    if u1 == "continuous":
        R = sy.Identity(2)

    if u1 == "const":
        R = sy.Identity(1)
        u = sy.Matrix([u[1]])         # take just u2

    if no_dead:
        Q = sy.Matrix([[0, 0, 0],
                       [0, .5, 0], 
                       [0, 0, 0]])
        
    else:
        Q=sy.Matrix([[0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 1]])
        
    return x.T@Q@x + u.T@R@u

def H(f, L, lam):
    """Define the hamiltonian for a system defined by dx = f, J[u] = int_0^inf L"""
    return lam @ f - L

def hamiltionian_partials(H, x, lam, u):
    ''' Takes in Sympy expressions for H, x, lam, 
    returns the state and costate eqns as well as Dh/Du
    as vectors of sympy expressions'''

    # make them each Sympy matrices
    x = sy.Matrix(x)
    lam = sy.Matrix(lam)

    # sympy can only do one derivative at a time, stack them into vectors
    x_dot = sy.Matrix([sy.diff(H, l) for l in lam])     
    lam_dot = -sy.Matrix([sy.diff(H, xi) for xi in x])

    # Can't remember if we'll use Dh/Du = 0
    zero = sy.Matrix([sy.diff(H, ui) for ui in u])

    return x_dot, lam_dot, zero 



### SYMPY MODEL DEFINITION
# define model
def sy_f_full(seasonal=True, no_dead=True):
    """This is the method to generate a full formulation of our model.

    - non-seasonal: autonomous, so we can use infinite horizon LQR
    - seasonal: non-autonomous, so infinite horizon LQR is no good
    
    Substitutes the dynamic variables into sy_f_simple. This can be used to actually
    evaluate when parameter values are known, but is much more painful to visualize than
    sy_f_simple."""
    # get dynamic parameter names and values
    names = sy_params_dynamic_names()
    tosub = sy_params_dynamic(seasonal=seasonal)

    # get simple f and sub in params
    if no_dead:
        fsimp = sy_f_simple_no_dead()
    else:
        fsimp = sy_f_simple()
    f = fsimp.subs(dict(zip(names, tosub)))

    return f

def sy_f_simple_no_dead():
    """Uses the sympy definitions below to construct a sympy expression for f(x)
    where x' = f(x) is our ODE. **This expression should only be used to construct
    a dynamic expression if any of the parameters are dynamic, and to double check the equations
    look right**. Can be used to construct both seasonal and autonomous versions of the ODE

    None of the sympy variables are functions in the output - and some of the parameters
    depend on t, u1, u2, and the state, so if you want to take derivatives or evaluate, sub in the dynamic
    expressions first.
    """
    # Removed the D class from the model
    S, I, T1, _ = sy_vars_model()
    alpha, beta, theta, xi, K = sy_params_dynamic_names()
    gamma, zeta, eta = sy_params_static()
    f = sy.Matrix([                 # f = [S', I', T1', D']^T
        [- alpha * I * S + gamma * T1 + eta * xi * I - zeta * S * (1 - T1 / K)],
        [alpha * I * S - zeta * I * (1 - T1 / K) - eta * xi * I - beta * I],
        [zeta * (I + S) * (1 - T1 / K) - gamma * T1]
    ])
    
    return f


### DEFINE CLASSES FOR DIFFERENT METHODS

class SolveBVP():
    """A class that uses solve_bvp to solve the axolotl problem"""

    def __init__(self, const_params, f, tf, initial_vals_x, initial_vals_lam=np.ones(3), no_dead=True, u1_const=150):
        """
        Parameters:
        - const_params (dict): dictionary of parameters (gamma, zeta, eta, alpha, beta)
        - f (function): sympy representation of the state evolution equation
        - tf (int): final time step
        - initial_vals (list): list of initial values for solve_bvp
        - no_dead (bool, opt): indicates whether to include the dead class
        - u1_const (float, opt): if present, indicates that the u1 control is constant and gives its value
        """
        self.const_params = const_params
        self.f = f
        self.tf = tf
        self.initial_x, self.initial_lam = initial_vals_x, initial_vals_lam
        self.no_dead = no_dead
        self.u1_const = u1_const

        # if u1_const, update constant parameters dictionary
        if not u1_const is None:
            self.const_params["u1"] = u1_const

    def _setup(self):
        """set up the problem
        
        Returns:
        - x_dot2: state transition equation with control substituted using PMP
        - lam_dot2: state transition equation with control substituted using PMP
        - u_symbolic: symbolic representation of u (used in substitutions)"""
        # define variables, lagrangian, and hamiltonian
        x, u, lam = define_variables()
        lagrangian = L(x, u)
        hamiltonian = H(self.f, lagrangian, lam)

        # get hamilton's equations
        if not self.u1_const is None:
            x_dot, lam_dot, stationary = hamiltionian_partials(hamiltonian, x, lam, sy.Matrix([u[1]]))
            u_sol = sy.simplify(sy.solve(stationary, u[1]))

            # u_sol is a dictionary?, so get the value
            u2 = sy.symbols("u2")
            u_symbolic = u_sol[u2]
        else:
            # this version needs some work but is probably the wrong model so we'll leave it for now
            x_dot, lam_dot, stationary = hamiltionian_partials(hamiltonian, x, lam, u)
            u_sol = sy.simplify(sy.solve(stationary, u))

        # plug in for u
        x_dot2 = sy.simplify(x_dot.subs(u_sol))
        lam_dot2 = sy.simplify(lam_dot.subs(u_sol))

        # store sympy variables as attributes
        self.x, self.u, self.lam = x, u, lam
        return x_dot2, lam_dot2, u_symbolic

    def _substitutions(self):
        """Convert sympy expressions to lambda functions
        
        Returns:
        - x_dot_func
        - lam_dot_func
        - u_func
        """
        x_dot, lam_dot, u = self._setup()

        # substitute in constant parameters
        x_dot_subs = x_dot.subs(self.const_params)
        lam_dot_subs = lam_dot.subs(self.const_params)
        u_subs = u.subs(self.const_params)

        # create lambda functions
        x_dot_func = sy.lambdify((self.x, self.lam), x_dot_subs, modules="numpy")
        lam_dot_func = sy.lambdify((self.x, self.lam), lam_dot_subs, modules="numpy")
        u_func = sy.lambdify((self.x, self.lam), u_subs, modules="numpy")

        return x_dot_func, lam_dot_func, u_func
    
    def solve(self):
        """Use solve_bvp to solve the system. Returns whatever solve_bvp returns."""
        x_dot_func, lam_dot_func, _ = self._substitutions()

        def ode(t, y):
            """Define an ode where the first n elements are the derivative the compartments
            and the last n elements are the derivative of the lambdas"""
            n = int(len(y)/2)                       # should work with or without the dead term
            dx = x_dot_func(y[:n], y[n:])
            dlam = lam_dot_func(y[:n], y[n:])

            return np.squeeze(np.concatenate((dx, dlam)))

        init_x = self.initial_x             
        init_lam = self.initial_lam 

        # define boundary conditions
        def bc(ya, yb):
            return np.squeeze(np.concatenate((ya[:len(init_x)] - init_x, yb[len(init_x):])))

        # get other necessary things
        tf = self.tf
        t_vals = np.linspace(0, tf, tf)
        final_x_guess = np.array([150, 200, 150])   # TODO: check this out to see if it is necessary/how much it changes things
        state_init = np.array([np.linspace(init_x[i], final_x_guess[i], len(t_vals)) for i in range(len(init_x))])
        costate_init = np.array([np.linspace(.1, .0001, len(init_x)) for _ in range(len(t_vals))]).T   # from Claude, suggested linear decay for costate initial conditions
        y_init = np.vstack((state_init, costate_init))

        return solve_bvp(ode, bc, t_vals, y_init)
    
    def plot_states(self, show=False):
        """Make a plot of the solution states"""
        sol = self.solve()
        plt.plot(sol.y[0], label="S")
        plt.plot(sol.y[1], label="I")
        plt.plot(sol.y[2], label="T1")
        plt.plot(500 - sol.y[0] - sol.y[1] - sol.y[2], label="D")
        plt.legend()
        plt.title("State Solutions")
        plt.xlabel("Day")
        plt.ylabel("Number of Frogs")

        if show:
            plt.show()

    def plot_u(self, show=False):
        u_function = self._substitutions()[-1]
        sol = self.solve()

        # plot the control parameters
        plt.plot(u_function(sol.y[:3], sol.y[3:]))
        plt.title("Fungal Bath Treatment Spending")
        plt.xlabel("Days")
        plt.ylabel("Amount Spent")

        if show:
            plt.show()


class SolveLQR():
    """A class that uses solve_bvp to solve the axolotl problem"""

    def __init__(self, const_params, f, no_dead=True, u1_const=0.1):
        """
        Parameters:
        - const_params (dict): dictionary of parameters (gamma, zeta, eta, alpha, beta)
        - f (function): sympy representation of the state evolution equation
        - no_dead (bool, opt): indicates whether to include the dead class
        - u1_const (float, opt): if present, indicates that the u1 control is constant and gives its value
        """
        self.const_params = const_params
        self.no_dead = no_dead
        self.f = f
        self.model = m.Model(seasonal=False, no_dead=True, const_u1=u1_const) # constant sauna maintenance

        # if u1_const, update constant parameters dictionary
        if u1_const:
            self.const_params["u1"] = u1_const


    def _solve(self, y0=np.array([475, 25, 0]), u0=np.array([1]), tf=5, t_steps=500):
        """set up the problem"""
        # define cost functional Matrixes
        self.Q = np.array([[0, 0, 0],   # cost of infected frog
                      [0, 3, 0],
                      [0, 0, 0]])

        self.R = np.eye(1)              # cost of control

        t_space = np.linspace(0, tf, t_steps)

        state_sol, control_sol, t_eval = self.iterated_lqr(lambda t, x, u: self.model.f_x(t, x, u), lambda t, x, u: self.model.f_u(t, x, u), t_space, y0, u0)

        return state_sol, control_sol, t_eval

    def solve_linearized_infinite(self, A, B, y0, t0, tf):
        '''
        Solve for control and state given a linearized system on an interval of length tn
        x_dot = Ax + Bu
        y0 - initial conditions (or last value or previous solution to linearized system)
        tn - length of one time that system is linearized on

        returns a solution to linerized system for tn years (run solve_linearized multiple times and then
        concatenate the solutions together)

        Basically one iteration of iterated LQR
        '''

        # solve continuous algebraic Riccati equation to get P
        P = la.solve_continuous_are(A, B, self.Q, self.R)


        # state equation Ax + 1/2 B R^-1 B^T lambda but lambda = -2Px by 
        def chytrid_evolve(t, y):
            return (A @ y - B @ B.T @ P @ y) # no R^-1 since thats just I

        # solve the optimal state evolution, and plot the position & velocity
        sol = solve_ivp(chytrid_evolve, [0, tf], y0, dense_output=True)
        t = np.linspace(t0, tf, 10)

        # S, I, T1, D = sol.sol(t)
        S = sol.sol(t)[0]
        I = sol.sol(t)[1]
        T1 = sol.sol(t)[2]

        u = - B.T @ P @ sol.sol(t)  # control

        return np.array([S[0], I[0], T1[0]]), u[:, 0]  # return first step of solution at point of linearization
    

    # this function iteratively calls solve_linearized_infinite
    def iterated_lqr(self, Df_x, Df_u, t_space, x0, u0):
        '''
        Solve for control and state solutions at discrete steps using LQR at each step
        Args:
            Df_x: a function to linearize with respect to x, returns matrix A for LQR
            Df_u: a function to linearize with respect to x, returns matrix B for LQR
            t_space: time steps on which to evaluate our solution
            x0: initial state
            u0: initial control

        Returns:
            state_solutions: array of optimal state over time
            control_solutions: array of optimal control over time
            acceptable_times: times where a finite P was found

        '''
        # store initial states in our solution
        state_solutions = [x0]
        control_solutions = [u0]
        acceptable_times = [t_space[0]]

        for i in range(len(t_space) - 1):   # iterate through each time step
            # get linearizations at this time step
            A = Df_x(0, state_solutions[-1], control_solutions[-1]) 
            B = Df_u(0, state_solutions[-1], control_solutions[-1])

            # uodate next step using LQR
            try:
                x_new, u_new = self.solve_linearized_infinite(A, B, state_solutions[-1], t_space[i], t_space[i+1])   # get next state using LQR
                state_solutions.append(x_new)
                control_solutions.append(u_new)
                acceptable_times.append(t_space[i])
            except:
                continue

        return np.array(state_solutions), np.array(control_solutions), acceptable_times   # optimal states and optimal controls 





#--------------------------------------------------------------
#        Init Parameters to get LQR plots in the jupyter notebook
#--------------------------------------------------------------

# # instantiate model object (to be used for linearization)
# model = m.Model(seasonal=False, no_dead=True, const_u1=0.1) # constant sauna maintenance

# # get optimal solutions
# x, u2 = iterated_lqr(lambda t, x, u: model.f_x(t, x, u), lambda t, x, u: model.f_u(t, x, u), t_space, y0, u0)
# S = x[:, 0]
# I = x[:, 1]
# T1 = x[:, 2]

# # plot
# plt.title("Optimal State for Linearized Chytrid")
# plt.plot(t_space, S, label="Susceptible")
# plt.plot(t_space, I, label="Infected")
# plt.plot(t_space, T1, label="In Treatment 1")
# plt.plot(t_space, 500 - S - I - T1, label="Death")  # plotting death total 500 - all other groups
# plt.legend()
# plt.show()

# # formatting
# plt.title("Optimal Control for Linearized Chytrid")
# plt.plot(t_space, u2, label=rf'$u_{2}$: Antifungal bath spending') 
# plt.legend()