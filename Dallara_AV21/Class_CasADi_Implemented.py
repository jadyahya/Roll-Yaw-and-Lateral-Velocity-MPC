import numpy as np
from scipy.linalg import block_diag
import sympy as sym
import casadi as cas
from casadi import *
import time
import csv

class Dallara_WMPC:
    def __init__(self) -> None:
        self.lat_vel = 0

    def set_parameters(self, ca, cs, reff, tf, tr, lf, lr, ms, hs,
                 m, ixx, ksf, ksr, ls, g, bsf,bsr, iz, kus_d, rmax, Ts, Q, R, alphamax, Qmax, deltamax, N):
        self.ca = ca
        self.cs = cs
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
        self.OldQ = Q
        self.R = R
        self.alphamax = alphamax
        self.Qmax = Qmax
        self.deltamax = deltamax
        self.N = N

        self.ks = (ksf+ksr)/2
        self.bs = (bsf+bsr)/2
    
    def set_state(self, xd, yd, speed, sideslipangle, wheel_speed, mux, muy, W0, phir, x0, current_acc, current_pos, current_ttl_index, current_pitch_rate):
        self.xd = xd
        self.yd = yd
        self.speed = speed
        self.sideslipangle = sideslipangle
        self.wheel_speed = wheel_speed
        self.mux = mux
        self.muy = muy
        self.W0 = W0
        self.phir = phir # That's it right?
        self.x0 = x0
        self.x0[1] = -self.lat_vel
        self.current_acc = current_acc
        self.current_pos = current_pos
        self.current_ttl_index = current_ttl_index
        self.current_pitch_rate = current_pitch_rate

    def load_estimator(self,u):
        faero_f = (309.6*(u/67.056)**2)/2.2 # Based on curve fir at rrh = 12.5 mm and frh = 10 mm
        faero_r = (668.7*(u/67.056)**2)/2.2 # Based on curve fir at rrh = 12.5 mm and frh = 10 mm
        fz1 = self.m*self.g/4 -self.ksf*self.ls/4*np.sin(self.x0[3])-self.bsf*self.ls/4*self.x0[4]*np.cos(self.x0[3])+faero_f/2
        fz2 = self.m * self.g / 4 +  self.ksf * self.ls / 4 * np.sin(self.x0[3]) + self.bsf * self.ls / 4 * self.x0[4] * np.cos(self.x0[3])+faero_f/2
        fz3 = self.m * self.g / 4 - self.ksr * self.ls / 4 * np.sin(self.x0[3]) - self.bsr * self.ls / 4 * self.x0[
            4] * np.cos(self.x0[3])+faero_r/2
        fz4 = self.m * self.g / 4 + self.ksr * self.ls / 4 * np.sin(self.x0[3]) + self.bsr * self.ls / 4 * self.x0[
            4] * np.cos(self.x0[3])+faero_r/2
        fz0 = np.array([fz1, fz2, fz3, fz4])
        return fz0
    
    def slip_angle_estimator(self,u):
        slip_angle_1 = self.W0[1] - (self.x0[1] + self.lf * self.x0[2]) / u * np.ones((2,))
        slip_angle_2 = -(self.x0[1] - self.lr * self.x0[2]) / u * np.ones((2,))
        slip_angle = np.hstack((slip_angle_1, slip_angle_2))
        return slip_angle

    # def lateral_load_estimator(self,ca1,ca2,ca3,ca4,slip_angle):
    #     fy1 = ca1 * (self.W0[1] - slip_angle[0])
    #     fy2 = ca2 * (self.W0[1] - slip_angle[1])
    #     fy3 = ca3 * (- slip_angle[2])
    #     fy4 = ca4 * (- slip_angle[3])
    #     fy0 = np.array([fy1, fy2, fy3, fy4]).reshape(4, )
    #     return fy0

    def slip_ratio_estimator(self,u):
        slip_ratio = np.zeros(4,)
        for i in range (4):
            if self.wheel_speed[i] == 0 or u < 10:
                slip_ratio[i] = 0.01
            else:
              slip_ratio[i] = (self.reff * self.wheel_speed[i] - u) / (self.reff * self.wheel_speed[i])
        return slip_ratio

    def speed_estimation(self):
        v = self.speed*np.sin(self.sideslipangle)
        self.x0[1] = v
        u = self.speed*np.cos(self.sideslipangle)
        return u,v


    def DugoffLinearizer(self, slip_angle, slip_ratio, ca, cs, mu, Fz):
        if slip_angle == 0:
            slip_angle = 0.01
      
        lamb_bar = mu * Fz * (1 + slip_ratio) / (2 * ((cs * slip_ratio) ** 2 + (ca * np.tan(slip_angle)) ** 2) ** 0.5)

        if slip_angle != 0 or slip_ratio != 0:
            if lamb_bar < 1:
                fbar = lamb_bar * (2 - lamb_bar)
                c_bar = (Fz * ca * mu * (np.tan(slip_angle) ** 2 + 1)) / (
                            2 * (ca ** 2 * np.tan(slip_angle) ** 2 + cs ** 2 * slip_ratio ** 2) ** (1 / 2)) - (
                                    Fz * ca ** 3 * mu * np.tan(slip_angle) ** 2 * (np.tan(slip_angle) ** 2 + 1)) / (
                                    2 * (ca ** 2 * np.tan(slip_angle) ** 2 + cs ** 2 * slip_ratio ** 2) ** (3 / 2))
                '''
                p1 = np.tan(slip_angle) ** 2 + 1
                p3 = ca**2*np.tan(slip_angle)**2+cs**2*slip_ratio**2
                p2 = Fz*mu*(slip_ratio+1)/(2*np.sqrt(p3))-2
                c_bar = ca**3*Fz*np.tan(slip_angle)**2*p2*p1/(2*p3)**3/2-ca*Fz*mu*p2*p1/2/np.sqrt(p3)+ca**3*Fz**2*mu**2*np.tan(slip_angle)**2*p1*(slip_ratio+1)/4/p2**2
                '''


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
    
    def Sx_Su_Sw_Sd_Sf(self, A, B, E, D, F):

        nX = np.size(A,0)
        nU = np.size(B,1)
        nE = np.size(E,1)
        
        Sx = np.eye(nX)
        A_tmp = A
        for i in range(self.N):
            Sx = np.vstack((Sx, A_tmp))
            A_tmp = A_tmp @ A
        
        SxB = Sx @ B
        Su = np.zeros((nX*(self.N+1),nU*self.N))
        for j in range(self.N):
            Su_tmp = np.vstack((np.zeros((nX, nU)), SxB[:-nX,:]))
            Su[:, 8*j:8*j+nU] = Su_tmp.reshape(Su_tmp.shape[0], nU)
            SxB = Su_tmp
        
        SxE = Sx @ E
        Sw = np.zeros((nX*(self.N+1),nE*self.N))
        for k in range(self.N):
            Sw_tmp = np.vstack((np.zeros((nX, nE)), SxE[:-nX,:]))
            Sw[:, 8*k:8*k+nE] = Sw_tmp.reshape(Sw_tmp.shape[0], nE)
            SxE = Sw_tmp
        
        Sd = np.zeros((nX,1)).reshape((nX,1))
        D_tmp = D.reshape(nX,1)
        for i in range(self.N):
            Sd = np.vstack((Sd, D_tmp))
            D_tmp = A @ D_tmp + D
            
        Sf = np.zeros((nX,1)).reshape((nX,1))
        F_tmp = F.reshape(nX,1)
        for i in range(self.N):
            Sf = np.vstack((Sf, F_tmp))
            F_tmp = A @ F_tmp + F

        return Sx, Su, Sw, Sd, Sf

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
        start_time = time.time()
        #u,v = self.speed_estimation()
        u = self.speed
        psi_desired = np.arctan(self.yd / self.xd)
        if u != 0.0:
            delta1 = psi_desired - self.x0[1] / u
            u = np.max((u, 10.0))
        else:
            delta1 = psi_desired
            u = np.max((u,10))

        fz0 = self.load_estimator(u)
        slip_angle = self.slip_angle_estimator(u)
        slip_ratio = self.slip_ratio_estimator(u)
        self.phir = self.x0[3]
        print(self.phir)
        print(self.x0)
        # Tire model Linearization
        fy1, ca1, a1 = self.DugoffLinearizer(slip_angle[0], slip_ratio[0], self.ca[0], self.cs[0], self.mux, fz0[0])
        fy2, ca2, a2 = self.DugoffLinearizer(slip_angle[1], slip_ratio[1], self.ca[1], self.cs[1], self.mux, fz0[1])
        fy3, ca3, a3 = self.DugoffLinearizer(slip_angle[2], slip_ratio[2], self.ca[2], self.cs[2], self.mux, fz0[2])
        fy4, ca4, a4 = self.DugoffLinearizer(slip_angle[3], slip_ratio[3], self.ca[3], self.cs[3], self.mux, fz0[3])
        B11 = np.array([[0, 0, 0, 0, 0], [0, -ca1 / u, -self.lf * ca1 / u, 0, 0]])
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
        
        # den = self.m**2*self.hs**2+self.ms*self.m*self.hs**2+self.ixx*self.m
        # Af = np.array([[0,1,0,0,0],
        #                [0,0,-u,(0.5*self.hs*self.ks*self.ls**2*self.ms-self.g*self.hs**2*self.m*self.ms)/den, self.bs*self.hs*self.ls**2*self.ms/den],
        #                [0,0,0,0,0],
        #                [0,0,0,0,1],
        #                [0,0,0,(self.g*self.hs*self.m**2-self.ks*0.5*self.ls**2*self.m)/den, -self.bs*self.ls**2*self.m/2/den]])
        # Bf = np.array([[0,0,0],
        #                [0,self.ixx/den+self.hs**2*self.m/den,0],
        #                [0,0,1/self.iz],
        #                [0,0,0],
        #                [0,self.m*self.hs/den,0]])
        # Cpr = np.array([[0],
        #                 [-self.g*self.ixx*self.m/den-self.g*self.hs**2*self.m**2/den],
        #                 [0],
        #                 [0],
        #                 [-self.g*self.hs*self.m**2/den]])

        den = self.m**2*self.hs**2+self.ms*self.m*self.hs**2+self.ixx*self.m
        Af = np.array([[0,1,0,0,0],
                       [0,0,-u,(0.5*self.hs*self.ks*self.ls**2*self.ms-self.g*self.hs**2*self.m*self.ms)/den, 0.5*self.bs*self.hs*self.ls**2*self.ms/den],
                       [0,0,0,0,0],
                       [0,0,0,0,1],
                       [0,0,0,(self.g*self.hs*self.m**2-self.ks*0.5*self.ls**2*self.m)/den, -self.bs*self.ls**2*self.m/2/den]])
        Bf = np.array([[0,0,0],
                       [0,self.ixx/den+self.hs**2*self.m/den,0],
                       [0,0,1/self.iz],
                       [0,0,0],
                       [0,self.m*self.hs/den,0]])
        Cpr = np.array([[0],
                        [2*self.g*self.hs**2*self.m*self.ms/den+self.g*self.hs**2*self.m**2/den+self.g*self.ixx*self.m/den],
                        [0],
                        [0],
                        [-self.g*self.hs*self.m**2/den]])
                        
        '''
        Testing with model from book by Pacejka
        den = self.m**2*self.hs**2-self.ms*self.m*self.hs**2+self.ixx*self.m
        Af = np.array([[0,1,0,0,0],
                       [0,0,-u,(0.5*self.hs*self.ks*self.ls**2*self.ms-self.g*self.hs**2*self.m*self.ms)/den, self.bs*self.hs*self.ls**2*self.ms/den],
                       [0,0,0,0,0],
                       [0,0,0,0,1],
                       [0,0,0,(self.g*self.hs*self.m**2-self.ks*0.5*self.ls**2*self.m)/den, -self.bs*self.ls**2*self.m/2/den]])
        Bf = np.array([[0,0,0],
                       [0,self.ixx/den+self.hs**2*self.m/den,0],
                       [0,0,1/self.iz],
                       [0,0,0],
                       [0,-self.m*self.hs/den,0]])
        Cpr = np.array([[0],
                        [self.g*self.ixx*self.m/den+self.g*self.hs**2*self.m**2/den],
                        [0],
                        [0],
                        [-self.g*self.hs*self.m**2/den]])
                        '''
        # X_dot = AX+EW+BU+D+FPhir
        # W = wheel torque, steer angle for each wheel (8x1)
        # States are y v r phi, phi_dot (ignoring tyre dynamics)
        # U = change in wheel torque, change in steering angle for each tire
        Tw = block_diag(Tw1,Tw2,Tw3,Tw4)
        Lw = block_diag(Lw1, Lw2, Lw3, Lw4)
        B1 = np.hstack((B11.T, B12.T, B13.T, B14.T)).T
        B2 = block_diag(B21, B22, B23, B24)
        D1 = np.hstack((D11,D12,D13,D14)).T
        A = Af+Bf@Lc@Lw@B1
        E = Bf@Lc@Lw@B2
        B = Bf@Lc@Lw@B2@Tw
        D = Bf@Lc@Lw@D1
        D = D.reshape(5,1)
        F = Cpr
        ###############################################################################
        ## CONTROL OBJECTIVES ##
        rmax = self.muy*self.g/u
        l= self.lf+self.lr
        rb = u/(l+self.kus_d*u**2)*delta1
        
        rd = np.sign(delta1)*np.min([np.abs(rb), rmax])
        if isinstance(rd, float):
            rd = rd
        else: 
            rd =rd.item()
        yd = self.yd.__float__()
        rd = rd.__float__()
        xd = np.array([-yd,0, rd, self.phir ,0]) #Negating yd is due to a difference in axes, the reasons of which are unclear
        #xd = np.array([0,0, rd, self.phir ,0])
        
        xd = xd.reshape(5,1)
        nX = 5
        nU = 8
        nW = 8
        nphir = 1
        # In the paper they use zero order hold, I will use euler discretization
        Ad = self.Ts*A + np.eye(np.size(A,0))
        Bd = self.Ts*B
        Ed = self.Ts*E
        Dd = self.Ts*D
        Fd = self.Ts*F
        #Formulating Problem:
        if np.abs(yd)>3.0:
            self.Q[0,0] = 2
            self.Q[2,2] = 2 # This is what I changed right now. If we increase it won't it make more effort? Yes after it gets into yd range it puts that much effort from what we just saw right. I see
        else:
        #    self.Q[0,0] = 10*(2.75-np.abs(yd))+150
        #    self.Q[2,2] = 400*(2.75-np.abs(yd))+800
             self.Q[0,0] = 15*(3.0-np.abs(yd))+25
             self.Q[2,2] = 160*(3.0-np.abs(yd))+150
        Sx, Su, Sw, Sd, Sf = self.Sx_Su_Sw_Sd_Sf(Ad, Bd, Ed, Dd, Fd) ## Check blocks 2 and 3 below for explanation
        # We need a weighting matrix Q. I would say emphasis on tracking lateral position, then heading, then lateral velocity, with little concern for roll charact.
        PN = self.Q
        Q_temp = np.kron(np.eye(self.N),self.Q)
        Qbar  = block_diag(Q_temp,PN)
        Rbar = 1500* np.eye(nU*self.N)

        Wbar = np.kron(np.ones((1,self.N)), self.W0).T #W0 recieved from vehicle, vector of inputs at time zero #### check
        Xdbar = np.kron(np.ones((1,self.N+1)), xd.T).T
        # Check block 4
        H = Su.T@Qbar@Su + Rbar
        F1 = Sx.T @Qbar @ Su
        F2 = Sw.T @ Qbar @ Su
        F3 = Sd.T@Qbar@Su
        F4 = -Qbar@Su
        F5 = Sf.T@Qbar@Su

        # Finally Cost Matricies:
        P = H
        x0 = self.x0.reshape(5,1)
        q = x0.T@F1+Wbar.T@F2 +F3+ self.phir*F5+ Xdbar.T@F4 #x0 from vehicle, xdes from path planner
        # Ax= np.array([[0,0,1,0,0],[0,0,-1,0,0],[0,-1/u,self.lr/u,0,0],[0,1/u,-self.lr/u,0,0]])
        # bx = np.array([rmax, rmax, self.alphamax,self.alphamax]).reshape(4,1)
        # Af = np.array([[0,0,1,0,0],[0,0,-1,0,0],[0,-1/u,self.lr/u,0,0],[0,1/u,-self.lr/u,0,0]])
        # bf = np.array([rmax, rmax, self.alphamax,self.alphamax]).reshape(4,1)

        Ax= np.array([[0,-1/u,self.lr/u,0,0],[0,1/u,-self.lr/u,0,0]])
        bx = np.array([self.alphamax,self.alphamax]).reshape(2,1)
        Af = np.array([[0,-1/u,self.lr/u,0,0],[0,1/u,-self.lr/u,0,0]])
        bf = np.array([ self.alphamax,self.alphamax]).reshape(2,1)


        # Ax= np.array([[0,-1/u,self.lr/u,0,0],[0,1/u,-self.lr/u,0,0]])
        # bx = np.array([np.inf,np.inf]).reshape(2,1)
        # Af = np.array([[0,-1/u,self.lr/u,0,0],[0,1/u,-self.lr/u,0,0]])
        # bf = np.array([np.inf,np.inf]).reshape(2,1)

        Au2 = np.array([[0,1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1]])
        Au3 = np.array([[0,1,0,-1,0,0,0,0],[0,-1,0,1,0,0,0,0]])
        Au = np.vstack((Au2,-Au2,Au3))

        bufip = self.deltamax*np.array([1,1,1,1])-np.array([self.W0[1],self.W0[3],self.W0[5],self.W0[7]])
        busip = self.deltamax*np.array([1,1,1,1])+np.array([self.W0[1],self.W0[3],self.W0[5],self.W0[7]])
        bufip = bufip.reshape(4,1)
        busip = busip.reshape(4,1)
        bu = np.vstack((bufip,busip,0.0001,0.0001))

        # Check Block 5
        G0_firstpart= block_diag(np.kron(np.eye(self.N),Au))
        temp = block_diag(np.kron(np.eye(self.N),Ax),Af)
        G0_sp = temp@Su
        W_contrib = temp@Sw@Wbar
        D_contrib = temp@Sd
        Phir_cont = temp@Sf*self.phir
        G0 = np.vstack((G0_firstpart,G0_sp))
        w0_2 = np.vstack((np.kron(np.ones((self.N,1)),bx), bf))
        w0_2 = w0_2-W_contrib-D_contrib-Phir_cont
        w0 = np.vstack((np.kron(np.ones((self.N,1)),bu), w0_2))
        e0_sp = -Ax
        Ax_tmp = -Ax@Ad
        Ax_tent = -Ad
        for i in range(self.N-1):
            e0_sp = np.vstack((e0_sp, Ax_tmp))
            Ax_tmp = Ax_tmp @ Ad
            Ax_tent = Ax_tent @Ad
        e0_sp = np.vstack((e0_sp, Af@Ax_tent))
        e0_fp = np.zeros((self.N*np.shape(bu)[0], 5)) ##2 length of Bu *n
        E0 = np.vstack((e0_fp,e0_sp))
        A = G0
        b = E0@x0+w0
        b = b.reshape(np.shape(b)[0],)

        # Optimizing with CasADI
        qp ={}
        H  = 2*cas.DM(H)
        g = cas.DM(q)

        uba = cas.DM(b.astype(float))
        A = cas.DM(A.astype(float))
        A = A.reshape((A.size1(),H.size1()))

        qp ={}
        qp['h']= H.sparsity()
        qp['a']=A.sparsity()
        try:
            # opts = {'printLevel': 'none'}
            S = cas.conic('S', 'osqp', qp)
            r = S(h=H, g=q, a=A, uba=uba)
            uOpt = r['x']
            JOpt = r['cost']
        except Exception as e:
            print(e)
            uOpt = np.zeros((80,1))
        x_opt = Sx@x0+Su@uOpt+Sw@Wbar+Sd+Sf*self.phir 
        a = self.W0[1]+float(uOpt[1])
        self.lat_vel = x_opt[6]
        end_time = time.time()
        elapsed_time = end_time - start_time
        a11 = [float(x_opt[5]),float(x_opt[6]),float(x_opt[7]),float(x_opt[8]),float(x_opt[9])]
        # with open('example.csv', mode='a', newline='') as file:
        #     # Create a writer object
        #     writer = csv.writer(file)
            
        #     # Write the header row
        #     writer.writerow(self.x0)
        #     writer.writerow(self.xd)
        #     writer.writerow(self.yd)
        #     writer.writerow(self.speed)
        #     writer.writerows(self.wheel_speed)
        #     writer.writerow(self.yd/self.xd)
        #     writer.writerow(self.mux)
        #     writer.writerow(self.muy)
        #     writer.writerow(self.phir)


        #     """
        #     self.xd = xd
        #     self.yd = yd
        #     self.speed = speed
        #     self.sideslipangle = sideslipangle
        #     self.wheel_speed = wheel_speed
        #     self.mux = mux
        #     self.muy = muy
        #     self.W0 = W0
        #     self.phir = phir # That's it right?
        #     """

        #     writer.writerow(a11)
        print("Time elapsed:", elapsed_time, "seconds")
        print(ca1)
        print(ca2)
        print(ca3)
        print(ca4)
        return a
