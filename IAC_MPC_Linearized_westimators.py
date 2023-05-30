import numpy as np
from scipy.linalg import block_diag
import sympy as sym
import pyomo.environ as pyo
from pyomo.opt import SolverStatus, TerminationCondition


class Dallara_WMPC:
    def __init__(self, xd, yd, speed, sideslipangle, wheel_speed, ca, cs, mux, muy, reff, tf, tr, lf, lr, ms, hs,
                 m, ixx, ksf, ksr, ls, g, bsf,bsr, iz, kus_d, rmax, Ts, Q, R, W0, phir, alphamax, Qmax, deltamax, x0, N):
        self.xd = xd
        self.yd = yd
        self.speed = speed
        self.sideslipangle = sideslipangle
        self.wheel_speed = wheel_speed
        self.ca = ca
        self.cs = cs
        self.mux = mux
        self.muy = muy
        self.reff = reff
        self.tf = tf
        self.tr = tr
        self.lf = lf
        self.lr = lr
        self.ms = ms
        self.hs = hs
        self.m = m
        self.ixx = ixx
        self.ksf = ksf
        self.ksr = ksr
        self.ls = ls
        self.g = g
        self.bsf = bsf
        self.bsr = bsr
        self.iz = iz
        self.kus_d = kus_d
        self.rmax = rmax
        self.Ts = Ts
        self.Q = Q
        self.R = R
        self.W0 = W0
        self.phir = phir
        self.alphamax = alphamax
        self.Qmax = Qmax
        self.deltamax = deltamax
        self.x0 = x0
        self.N = N
        return

    def load_estimator(self):
        fz1 = self.m*self.g/4 +self.ksf*self.ls/4*np.sin(self.x0[3])+self.bsf*self.ls/4*self.x0[4]*np.cos(self.x0[3])
        fz2 = self.m * self.g / 4 -  self.ksf * self.ls / 4 * np.sin(self.x0[3]) - self.bsf * self.ls / 4 * self.x0[4] * np.cos(self.x0[3])
        fz3 = self.m * self.g / 4 + self.ksr * self.ls / 4 * np.sin(self.x0[3]) + self.bsr * self.ls / 4 * self.x0[
            4] * np.cos(self.x0[3])
        fz4 = self.m * self.g / 4 - self.ksr * self.ls / 4 * np.sin(self.x0[3]) - self.bsr * self.ls / 4 * self.x0[
            4] * np.cos(self.x0[3])
        fz0 = np.array([fz1, fz2, fz3, fz4])
        return fz0

    def slip_angle_estimator(self,u):
        slip_angle_1 = self.W0[1] - (self.x0[1] + self.lf * self.x0[2]) / u * np.ones((2,))
        slip_angle_2 = (self.x0[1] - self.lr * x0[2]) / u * np.ones((2,))
        slip_angle = np.hstack((slip_angle_1, slip_angle_2))
        return slip_angle

    def lateral_load_estimator(self,ca1,ca2,ca3,ca4,slip_angle):
        fy1 = ca1 * (self.W0[1] - slip_angle[0])
        fy2 = ca2 * (self.W0[1] - slip_angle[1])
        fy3 = ca3 * (- slip_angle[2])
        fy4 = ca4 * (- slip_angle[3])
        fy0 = np.array([fy1, fy2, fy3, fy4]).reshape(4, )
        return fy0

    def slip_ratio_estimator(self,u):
        for i in range(4):
            slip_ratio[i] = (self.reff * self.wheel_speed[i] - u) / np.max(self.reff * self.wheel_speed[i],u)
        return slip_ratio

    def speed_estimation(self):
        v = self.speed*np.sin(self.sideslipangle)
        self.x0[1] = v
        u = self.speed*np.sin(self.sideslipangle)
        return u,v


    def DugoffLinearizer(slip_angle, slip_ratio, ca, cs, mu, Fz):
        lamb_bar = mu * Fz * (1 + slip_ratio) / (2 * ((cs * slip_ratio) ** 2 + (ca * np.tan(slip_angle)) ** 2) ** 0.5)

        if slip_angle != 0 or slip_ratio != 0:
            if lamb_bar < 1:
                fbar = lamb_bar * (2 - lamb_bar)
                c_bar = (Fz * ca * mu * (np.tan(slip_angle) ** 2 + 1)) / (
                            2 * (ca ** 2 * np.tan(slip_angle) ** 2 + cs ** 2 * slip_ratio ** 2) ** (1 / 2)) - (
                                    Fz * ca ** 3 * mu * np.tan(slip_angle) ** 2 * (np.tan(slip_angle) ** 2 + 1)) / (
                                    2 * (ca ** 2 * np.tan(slip_angle) ** 2 + cs ** 2 * slip_ratio ** 2) ** (3 / 2))


            else:
                fbar = 1
                c_bar = (ca * (np.tan(slip_angle) ** 2 + 1)) / (slip_ratio + 1)
            Fy_bar = ca * np.tan(slip_angle) * fbar / (1 + slip_ratio)
            c_bar = c_bar.item()
            Fy_bar = Fy_bar.item()

        else:
            fbar = 1
            c_bar = ca.item()
            Fy_bar = 0

        return Fy_bar, c_bar, slip_angle

    def ctrl_action(self):
        # States: y,v,psi, r, phi, phi_dot
        # Inputs: Change in each wheel torque and in each wheel steer angle
        """
        STATES:
        y: lateral position
        vy: lateral speed
        psi: heading angle
        r: yaw rate
        phi: roll angle
        phid: roll rate

        =================================================================
        =================================================================
        Input Definitions:

        slip_angle: tire slip angle at each tire (aka alpha) --- (5x1)
        slip_ratio: tire slip ratio at each tire (aka sigma) --- (5x1)
        ca: tire cornering stifness (aka C_alpha)            --- (5x1)
        cs: tire longitudinal stiffness (aka C_sigma)        --- (5x1)
        mux: coefficient of friction in x                    --- (scalar)
        muy: coefficient of friction in y                    --- (scalar)
        fz0: tire vertical load                              --- (5x1)
        fyo: tire lateral load                               --- (5x1)
        u: longitudinal velcoity                             --- (scalar)
        reff: tire effective radius of rotation              --- (scalar)
        delta1: initial commanded steering input             --- (scalar)
        tf: front wheel track (dist between tire centers)    --- (scalar)
        tr: rear wheel track (dist between tire centers)     --- (scalar)
        lf: distance from cg to front wheels                 --- (scalar)
        lr: distance from cg to rear wheels                  --- (scalar)
        ms: vehicle sprung mass (without mass of tires and suspension) --- (scalar)
        hs: cg of sprung mass to roll center                 --- (scalar)
        m: total mass of vehicle (sprung and unsprung)       --- (scalar)
        ixx: moment of inertia about the x-axis (roll)       --- (scalar)
        ks: spring stiffness of suspension                   --- (scalar)
        ls: distance between suspension fixtures to vehicle  --- (scalar)
        g: gravity                                           --- (scalar)
        bs: suspension damping coefficient                   --- (scalar)
        iz: vehivle moment of inertia about z axis (yaw)     --- (scalar)
        kus_d: desired understeer coefficient of vehicle     --- (scalar)
        rmax: maximum allowed yaw rate (to be selected)      --- (scalar)
        Ts: sampling time for euler discretization           --- (scalar)
        Q: cost matrix for states                            --- (5x5)
        R: cost matrix for inputs                            --- (5x5)
        W0: initial torques and steering angles [Q1,delta1,Q2,delta2,Q3,delta3,Q4,delta4] ---(8x1)
        Phir: banking angle at over the prediction horizon   --- (Nx1)
        alphamax: maximum allowable rear wheel side slip     --- (scalar)
        Qmax: maximum torque attainable at each wheel        --- (scalar)
        deltamax: maximum steering angle                     --- (scalar)
        x0: initial states [y,vy,r,phi,phid]                 --- (5x1)
        psid: desired heading angle                          --- (scalar)
        yd:desired lateral position                          --- (scalar)
        N: MPC Prediction horizon                            --- (scalar)
        """

        u,v = speed_estimation(self)
        psi_desired = np.arctan(self.yd / self.xd)
        if u != 0.0:
            delta1 = psi_desired - self.x0[1] / u
        else:
            delta1 = psi_desired

        fz0 = load_estimator(self)
        slip_angle = slip_angle_estimator(self)
        slip_ratio = slip_ratio_estimator(self)

        # Tire model Linearization
        fy1, ca1, a1 = self.DugoffLinearizer(self.slip_angle[0], self.slip_ratio[0], self.ca[0], self.cs[0], self.mux, self.fz0[0])
        fy2, ca2, a2 = self.DugoffLinearizer(self.slip_angle[1], self.slip_ratio[1], self.ca[1], self.cs[1], self.mux, self.fz0[1])
        fy3, ca3, a3 = self.DugoffLinearizer(self.slip_angle[2], self.slip_ratio[2], self.ca[2], self.cs[2], self.mux, self.fz0[2])
        fy4, ca4, a4 = self.DugoffLinearizer(self.slip_angle[3], self.slip_ratio[3], self.ca[3], self.cs[3], self.mux, self.fz0[3])
        B11 = np.array([[0, 0, 0, 0, 0], [0, -self.ca1 / u, -self.lf * ca1 / u, 0, 0]])
        B21 = np.array([[1 / self.reff, 0], [0, ca1]])
        D11 = np.array([0, fy1 - ca1 * a1])

        B12 = np.array([[0, 0, 0, 0, 0], [0, -ca2 / u, -self.lf * ca2 / u, 0, 0]])
        B22 = np.array([[1 / self.reff, 0], [0, ca2]])
        D12 = np.array([0, fy2 - ca2 * a2])

        B13 = np.array([[0, 0, 0, 0, 0], [0, -ca3 / u, self.lr * ca3 / u, 0, 0]])
        B23 = np.array([[1 / self.reff, 0], [0, ca3]])
        D13 = np.array([0, fy3 - ca3 * a3])

        B14 = np.array([[0, 0, 0, 0, 0], [0, -ca4 / u, self.lr * ca4 /u, 0, 0]])
        B24 = np.array([[1 / self.reff, 0], [0, ca4]])
        D14 = np.array([0, fy4 - ca4 * a4])

        # Local actuator reconfiguration matrix: 0 inactive actuator 1 active actuator
        tq1 = 0
        tq2 = 0
        tq3 = 0
        tq4 = 0
        tdelta1 = 1
        tdelta2 = 1
        tdelta3 = 0
        tdelta4 = 0

        Tw1 = np.array([[tq1, 0], [0, tdelta1]])
        Tw2 = np.array([[tq2, 0], [0, tdelta2]])
        Tw3 = np.array([[tq3, 0], [0, tdelta3]])
        Tw4 = np.array([[tq4, 0], [0, tdelta4]])

        # Reconfiguration (Mapping) Matricies: Delta1 input
        Lw1 = np.array([[np.cos(self.W0[1]), -np.sin(self.W0[1])], [np.sin(self.W0[1]), np.cos(self.W0[1])]])
        Lw2 = np.array([[np.cos(self.W0[1]), -np.sin(self.W0[1])], [np.sin(self.W0[1]), np.cos(self.W0[1])]])
        Lw3 = np.array([[1, 0], [0, 1]])
        Lw4 = np.array([[1, 0], [0, 1]])

        # Mapping Matrix from corner forces
        Lc = np.array([[1, 0, 1, 0, 1, 0, 1, 0], [0, 1, 0, 1, 0, 1, 0, 1],
                       [-self.tf / 2, self.lf, self.tf / 2, self.lf, -self.tr / 2, -self.lr, self.tr / 2, -self.lr]])

        # Vehicle Body Dynamics y vy psi,psidot roll, rollrate
        # Equations written in block 1 below
        den = 2 * self.m * (self.ixx + self.hs ** 2 * self.m + self.hs ** 2 * self.ms)
        Af = np.array([[0, 1, 0, 0, 0],
                       [0, 0, -u, (self.hs * self.ks * self.ls ** 2 * self.ms + 2 * self.g * self.ixx * self.m + 2 * self.g * self.hs ** 2 * self.m ** 2) / den,
                        self.bs * self.hs * self.ls ** 2 * self.ms / den], [0, 0, 0, 0, 0], [0, 0, 0, 0, 1],
                       [0, 0, 0, (2 * self.g * self.hs * self.m * self.m - self.ks * self.ls ** 2 * self.m) / den, -self.bs * self.ls * self.ls * self.m / den]])
        Bf = np.array(
            [[0, 0, 0], [0, 2 * (self.ixx + self.m * self.hs ** 2) / den, 0], [0, 0, 1 / self.iz], [0, 0, 0], [0, 2 * self.m * self.hs / den, 0]])
        Cpr = np.array([[0], [(2 * self.g * self.hs ** 2 * self.m * self.ms - self.hs * self.ks * self.ls ** 2 * self.ms) / den], [0], [0],
                        [(self.ks * self.ls ** 2 * self.m - 2 * self.g * self.hs * self.m * self.m) / den]])

        # X_dot = AX+EW+BU+D+FPhir
        # W = wheel torque, steer angle for each wheel (8x1)
        # States are y v r phi, phi_dot (ignoring tyre dynamics)
        # U = change in wheel torque, change in steering angle for each tire

        Tw = block_diag(Tw1, Tw2, Tw3, Tw4)
        Lw = block_diag(Lw1, Lw2, Lw3, Lw4)
        B1 = np.hstack((B11.T, B12.T, B13.T, B14.T)).T
        B2 = block_diag(B21, B22, B23, B24)
        D1 = np.hstack((D11, D12, D13, D14)).T
        A = Af + Bf @ Lc @ Lw @ B1
        E = Bf @ Lc @ Lw @ B2
        B = Bf @ Lc @ Lw @ B2 @ Tw
        D = Bf @ Lc @ Lw @ D1
        D = D.reshape(5, )
        F = Cpr

        ###############################################################################
        ## CONTROL OBJECTIVES ##
        rmax = self.muy * self.g / u
        l = self.lf + self.lr
        rb = u / (self.l + self.kus_d * u ** 2) * delta1
        rd = np.sign(delta1) * np.min([np.abs(rb), rmax])
        if isinstance(rd, float):
            rd = rd
        else:
            rd = rd.item()
        xd = np.array([self.yd, 0, rd, self.phir, 0])
        xd = xd.reshape(5, )
        nX = 5
        nU = 8
        nW = 8
        nphir = 1
        # In the paper they use zero order hold, I will use euler discretization
        Ad = self.Ts * A + np.eye(np.size(A, 0))
        Bd = self.Ts * B
        Ed = self.Ts * E
        Dd = self.Ts * D
        Fd = self.Ts * F
        model = pyo.ConcreteModel()
        model.tidx = pyo.Set(initialize=range(0, N + 1))
        model.tidu = pyo.Set(initialize=range(0, N))
        model.xidx = pyo.Set(initialize=range(0, nX))
        model.uidx = pyo.Set(initialize=range(0, nU))
        model.x = pyo.Var(model.xidx, model.tidx)
        model.u = pyo.Var(model.uidx, model.tidx)
        model.phir = self.phir
        model.xd = xd
        model.Q = Q
        model.R = R
        model.N = self.N
        model.A = Ad
        model.B = Bd
        model.E = Ed
        model.D = Dd
        model.F = Fd
        model.W0 = W0
        model.P = Q

        # Objective:

        def objective_rule(model):
            costX = 0.0
            costU = 0.0
            costT = 0.0
            for t in model.tidx:
                for i in model.xidx:
                    for j in model.xidx:
                        if t < N:
                            costX += (model.x[i, t] - xd[i].item()) * model.Q[i, j] * (model.x[j, t] - xd[j].item())
            for t in model.tidx:
                for i in model.uidx:
                    for j in model.uidx:
                        if t < N:
                            costU += model.u[i, t] * model.R[i, j] * model.u[j, t]
            for i in model.xidx:
                for j in model.xidx:
                    costT += model.x[i, model.N] * model.P[i, j] * model.x[j, model.N]
            return costX + costU + costT

        model.Obj = pyo.Objective(rule=objective_rule, sense=pyo.minimize)
        '''
        print(xd[0])
        print(xd[1])
        print(xd[2])
        print(xd[3])
        print(xd[4])
        model.Obj = pyo.Objective(expr = (sum(sum(model.R[j,j]*(model.u[j, t])**2 for j in model.uidx) for t in model.tidu ) + sum(sum(model.Q[i,i]*(model.x[i, t]-xd[i].item())**2 for i in model.xidx) for t in model.tidx)) , sense=pyo.minimize)
        '''

        # Model Dynanamics
        def eq_const_rule(model, i, t):
            return model.x[i, t + 1] - (sum(model.A[i, j] * model.x[j, t] for j in model.xidx) + sum(
                model.E[i, j] * model.W0[j] for j in model.uidx) + model.phir * model.F[i].item() + sum(
                model.B[i, j] * model.u[j, t] for j in model.uidx) + model.D[
                                            i]) == 0 if t < model.N else pyo.Constraint.Skip

        model.constraint1 = pyo.Constraint(model.xidx, model.tidx, rule=eq_const_rule)

        # Input Constraints:

        '''
        TO BE ADDED IF WE OPT FOR LONG AND LAT CONTROL (WITH SOME MODIFICATIONS)
        # Torque Saturation:
        model.constraint2 = pyo.Constraint(model.tidx, rule = lambda model, t: model.u[0,t] <= Qmax - model.W0[0] if t < N else pyo.Constraint.Skip)
        model.constraint3 = pyo.Constraint(model.tidx, rule = lambda model, t: model.u[2,t] <= Qmax - model.W0[2] if t < N else pyo.Constraint.Skip)
        model.constraint4 = pyo.Constraint(model.tidx, rule = lambda model, t: model.u[4,t] <= Qmax - model.W0[4] if t < N else pyo.Constraint.Skip)
        model.constraint5 = pyo.Constraint(model.tidx, rule = lambda model, t: model.u[6,t] <= Qmax - model.W0[6] if t < N else pyo.Constraint.Skip)

        model.constraint6 = pyo.Constraint(model.tidx, rule = lambda model, t: model.u[0,t] >= Qmin - model.W0[0] if t < N else pyo.Constraint.Skip)
        model.constraint7 = pyo.Constraint(model.tidx, rule = lambda model, t: model.u[2,t] >= Qmin - model.W0[2] if t < N else pyo.Constraint.Skip)
        model.constraint8 = pyo.Constraint(model.tidx, rule = lambda model, t: model.u[4,t] >= Qmin - model.W0[4] if t < N else pyo.Constraint.Skip)
        model.constraint9 = pyo.Constraint(model.tidx, rule = lambda model, t: model.u[6,t] >= Qmin - model.W0[6] if t < N else pyo.Constraint.Skip)


        fx0_p=np.zeros((4,1))
        for i in range(4):
          fx0_p[i] = mux*fz0[i]*np.sqrt(1-(fy0[i]/(muy*fz0[i]))**2)

        # Friction Circle
        model.constraint10 = pyo.Constraint(model.tidx, rule = lambda model, t: model.u[0,t] <= fx0_p[0].item() - model.W0[0] if t < N else pyo.Constraint.Skip)
        model.constraint11= pyo.Constraint(model.tidx, rule = lambda model, t: model.u[2,t] <= fx0_p[1].item() - model.W0[2] if t < N else pyo.Constraint.Skip)
        model.constraint12 = pyo.Constraint(model.tidx, rule = lambda model, t: model.u[4,t] <= fx0_p[2].item() - model.W0[4] if t < N else pyo.Constraint.Skip)
        model.constraint13 = pyo.Constraint(model.tidx, rule = lambda model, t: model.u[6,t] <= fx0_p[3].item() - model.W0[6] if t < N else pyo.Constraint.Skip)

        model.constraint14 = pyo.Constraint(model.tidx, rule = lambda model, t: model.u[0,t] >= -fx0_p[0].item() - model.W0[0] if t < N else pyo.Constraint.Skip)
        model.constraint15 = pyo.Constraint(model.tidx, rule = lambda model, t: model.u[2,t] >= -fx0_p[1].item() - model.W0[2] if t < N else pyo.Constraint.Skip)
        model.constraint16 = pyo.Constraint(model.tidx, rule = lambda model, t: model.u[4,t] >= -fx0_p[2].item() - model.W0[4] if t < N else pyo.Constraint.Skip)
        model.constraint17 = pyo.Constraint(model.tidx, rule = lambda model, t: model.u[6,t] >= -fx0_p[3].item() - model.W0[6] if t < N else pyo.Constraint.Skip)
        '''
        # Maximum steering angle constraints
        model.constraint18 = pyo.Constraint(model.tidx, rule=lambda model, t: model.u[1, t] <= self.deltamax - model.W0[
            1] if t < self.N else pyo.Constraint.Skip)
        model.constraint19 = pyo.Constraint(model.tidx, rule=lambda model, t: model.u[3, t] <= self.deltamax - model.W0[
            3] if t < self.N else pyo.Constraint.Skip)
        model.constraint20 = pyo.Constraint(model.tidx, rule=lambda model, t: model.u[5, t] <= self.deltamax - model.W0[
            5] if t < self.N else pyo.Constraint.Skip)
        model.constraint21 = pyo.Constraint(model.tidx, rule=lambda model, t: model.u[7, t] <= self.deltamax - model.W0[
            7] if t < self.N else pyo.Constraint.Skip)

        model.constraint22 = pyo.Constraint(model.tidx, rule=lambda model, t: model.u[1, t] >= -self.deltamax - model.W0[
            1] if t < self.N else pyo.Constraint.Skip)
        model.constraint23 = pyo.Constraint(model.tidx, rule=lambda model, t: model.u[3, t] >= -self.deltamax - model.W0[
            3] if t < self.N else pyo.Constraint.Skip)
        model.constraint24 = pyo.Constraint(model.tidx, rule=lambda model, t: model.u[5, t] >= -self.deltamax - model.W0[
            5] if t < self.N else pyo.Constraint.Skip)
        model.constraint25 = pyo.Constraint(model.tidx, rule=lambda model, t: model.u[7, t] >= -self.deltamax - model.W0[
            7] if t < self.N else pyo.Constraint.Skip)

        ## State Constraints
        model.constraint26 = pyo.Constraint(model.tidx, rule=lambda model, t: model.x[
                                                                                  2, t] <= rmax if t <= self.N else pyo.Constraint.Skip)
        model.constraint27 = pyo.Constraint(model.tidx, rule=lambda model, t: model.x[
                                                                                  2, t] >= -rmax if t <=self. N else pyo.Constraint.Skip)
        model.constraint28 = pyo.Constraint(model.tidx, rule=lambda model, t: -model.x[1, t] / u + self.lr / u * model.x[
            2, t] <= self.alphamax if t <= self.N else pyo.Constraint.Skip)
        model.constraint29 = pyo.Constraint(model.tidx, rule=lambda model, t: -model.x[1, t] / u + self.lr / u * model.x[
            2, t] >= -self.alphamax if t <= self.N else pyo.Constraint.Skip)

        model.constraint32 = pyo.Constraint(model.xidx, rule=lambda model, i: model.x[i, 0] == self.x0[i])
        model.constraint30 = pyo.Constraint(model.tidx, rule=lambda model, t: model.u[1, t] == model.u[
            3, t] if t < self.N else pyo.Constraint.Skip)

        # model.constraint31 = pyo.Constraint(model.xidx, rule = lambda model, i: model.x[0,N] == xd[0])

        results = pyo.SolverFactory('ipopt').solve(model).write()

        # Plotting
        # plot results
        #x1 = [pyo.value(model.x[0, 0])]
        #x2 = [pyo.value(model.x[1, 0])]
        #x3 = [pyo.value(model.x[2, 0])]
        #x4 = [pyo.value(model.x[3, 0])]
        #x5 = [pyo.value(model.x[4, 0])]

        #u1 = [pyo.value(model.u[0, 0])]
        u2 = [pyo.value(model.u[1, 0])]
        #u3 = [pyo.value(model.u[2, 0])]
        #u4 = [pyo.value(model.u[3, 0])]
        #u5 = [pyo.value(model.u[4, 0])]
        #u6 = [pyo.value(model.u[5, 0])]
        #u7 = [pyo.value(model.u[6, 0])]
        #u8 = [pyo.value(model.u[7, 0])]

        return u2

