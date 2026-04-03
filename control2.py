import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
from scipy import linalg as la
import model as m


def solve_linearized_infinite(A, B, y0, t0, tf):
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
    P = la.solve_continuous_are(A, B, Q, R)


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
def iterated_lqr(Df_x, Df_u, t_space, x0, u0):
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

    '''
    # store initial states in our solution
    state_solutions = [x0]
    control_solutions = [u0]

    for i in range(len(t_space) - 1):   # iterate through each time step
        # get linearizations at this time step
        A = Df_x(0, state_solutions[i], control_solutions[i]) 
        B = Df_u(0, state_solutions[i], control_solutions[i])

        # uodate next step using LQR
        x_new, u_new = solve_linearized_infinite(A, B, state_solutions[i], t_space[i], t_space[i+1])   # get next state using LQR
        state_solutions.append(x_new)
        control_solutions.append(u_new)

    return np.array(state_solutions), np.array(control_solutions)   # optimal states and optimal controls 


#--------------------------------------------------------------
#        Init Parameters to get plots in the jupyter notebook
#--------------------------------------------------------------

Q = np.array([[0, 0, 0],   # cost of infected frog
             [0, 3, 0],
             [0, 0, 0]])

R = np.eye(1)              # cost of control

# initial conditions
y0 = np.array([475, 25, 0])
u0 = np.array([1])  # init guess of optimal control on bath spending

# establish time space
tf = 5
t_steps = 500
t_space = np.linspace(0, tf, t_steps)

# instantiate model object (to be used for linearization)
model = m.Model(seasonal=False, no_dead=True, const_u1=0.1) # constant sauna maintenance

# get optimal solutions
x, u2 = iterated_lqr(lambda t, x, u: model.f_x(t, x, u), lambda t, x, u: model.f_u(t, x, u), t_space, y0, u0)
S = x[:, 0]
I = x[:, 1]
T1 = x[:, 2]

# plot
plt.title("Optimal State for Linearized Chytrid")
plt.plot(t_space, S, label="Susceptible")
plt.plot(t_space, I, label="Infected")
plt.plot(t_space, T1, label="In Treatment 1")
plt.plot(t_space, 500 - S - I - T1, label="Death")  # plotting death total 500 - all other groups
plt.legend()
plt.show()

# formatting
plt.title("Optimal Control for Linearized Chytrid")
plt.plot(t_space, u2, label=rf'$u_{2}$: Antifungal bath spending') 
plt.legend()