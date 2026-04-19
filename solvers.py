import numpy as np
from matplotlib import pyplot as plt
from scipy import linalg as la
from scipy.integrate import solve_ivp, solve_bvp
import model

class Solver():
    def __init__(self, const_u1):
        self.const_u1 = const_u1
        
        # define a model and timescale
        PARAMS = {
            "α": 0.000035,
            "β": 0.003
        }

        self.INITIAL_X = np.array([475, 25, 0])

        self.INITIAL_u2guess = np.array([50])

        self.Qmult = 1
        self.Rmult = 1

        finite = False

        self.tf = 1000
        self.tsteps = 10001
        self.tspan = np.linspace(0, self.tf, self.tsteps)

        self.mod = model.Model(
            seasonal=False,
            params=PARAMS,
            const_u1=const_u1,
            no_dead=True
        )

        self.Qbase = np.array([[0, 0, 0],
                    [0, 1, 0],
                    [0, 0, 0]])    # cost of frog death

        self.Rbase = np.eye(1)   # cost of control


    def get_u_infinite(self, A, B, x, Q, R):
        '''
        Solve for control given a linearized system
        x_dot = Ax + Bu
        x - state conditions (or last value or previous solution to linearized system)
        '''

        # solve continuous algebraic Riccati equation to get P
        P = la.solve_continuous_are(A, B, Q, R)

        u = - np.linalg.inv(R) @ B.T @ P @ x  # control

        return u

    def get_u_finite(self, A, B, x, Q, R, tcurr, tfinal, m=0):
        '''
        Solve for control given a linearized system
        x_dot = Ax + Bu
        x - state conditions (or last value or previous solution to linearized system)
        '''
        t_steps = 100
        t_space = np.linspace(tcurr, self.tf, t_steps)

        def P_evolve(t, p): # riccatti ode
            P = p.reshape(A.shape[0], A.shape[0])
            # print("P", P, "\nA", A, "\nB", B, "\nQ", Q, "\nR", R)
            return (-P @ A - A.T @ P - Q + P @ B @ np.linalg.inv(R) @ B.T @ P).flatten()
        
        pf = np.diag([0, m, 0]).flatten()
        #pf = np.zeros(A.shape[0]**2) # P(t_f) = M = 0
        sol_P = solve_ivp(P_evolve, [tfinal, tcurr], pf, dense_output=True, t_eval=t_space[::-1])
        p = sol_P.y
        P = p.reshape(sol_P.t.shape[0], A.shape[0], A.shape[0])[::-1]

        u = - np.linalg.inv(R) @ B.T @ P[0] @ x  # solve for control

        return u

    def combinestate(self, x, u):
        return np.concatenate(x, u)

    def breakstate(self, y):
        return y[:3], y[-1]

    def control(self, t, x, uprev, Q, R, tfinal, finite=False):
        "Get the appropriate control using a previous guess"
        A = self.mod.f_x(t, x, uprev)
        B = self.mod.f_u(t, x, uprev)
        if not finite:
            return self.get_u_infinite(A, B, x, Q, R)
        else:
            return self.get_u_finite(A, B, x, Q, R, t, tfinal)

    def evolvex(self, t, x, u):
        "Evolve the system"
        # print("Evolving")
        # print(t, x, u)
        # print(mod.f(t, x, u))
        return self.mod.f(t, x, u)

    def solvelqr(self, uinit, xinit, Q, R, tf, finite=False, uupdateshrink=0.5):
        """Solve LQR"""
        constusteps = 10
        tspan_steps = [self.tspan[tidx:min(self.tsteps-1, tidx+constusteps)+1] for tidx in np.arange(0, self.tsteps-1, constusteps)]

        u = [uinit]
        x = self.INITIAL_X

        tvals = list()
        xvals = list()
        uvals = list()

        for steps in tspan_steps:
            t = steps[0]
            # update u based on current state and previous u value
            
            unewprop = self.control(t, x, u, Q, R, tf, finite=finite)
            u = (uupdateshrink * (unewprop-u)) + u
            # u = control(t, x, u, Q, R, tf, finite=finite)
            # possibly use newton's method?

            # update t, u, x
            tvals.append(t)
            xvals.append(x)
            uvals.append(u)

            # print(evolvex(
            #     t,
            #     x,
            #     u,
            # ))
            
            # solve on each interval, fixing u
            tspn = (steps[0], steps[-1])
            sol = solve_ivp(
                self.evolvex,
                tspn,
                x,
                args=(u,),
                t_eval=np.linspace(*tspn, 100)
            )

            x = sol.y[:, -1]

        return np.array(tvals), np.stack(xvals), np.array(uvals) 

    def costben(self, uinit, xinit, Q, R, tf, finite=False, uupdateshrink=0.5):
        ts, xs, us = self.solvelqr(uinit, self.INITIAL_X, Q, R, tf, finite=finite, uupdateshrink=uupdateshrink)
        # plt.plot(ts, xs)
        # plt.show()
        # plt.plot(ts, us)
        # plt.show()
        return xs[-1], (self.const_u1 * tf) + (np.sum(us) * len(us) / tf)

    def costbensweep(self):
        finalxs = list()
        totalcosts = list()

        for Rmult in np.logspace(-3, 3, 15, base=10):
            finalx, totalcost = self.costben(100, self.INITIAL_X, self.Qmult * self.Qbase, Rmult * self.Rbase, self.tf, finite=False)
            finalxs.append(finalx)
            totalcosts.append(totalcost)

        xvals = np.stack(finalxs)
        totalcosts = np.array(totalcosts)
        survived = np.sum(xvals, axis=1)

        return totalcosts, survived
