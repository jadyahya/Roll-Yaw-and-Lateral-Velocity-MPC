'''
delta: steering angle
vy : Lateral velocity
vx : longitudinal velocity
lf: distance from cg to front axle
r: yaw rate
lr: distance from cg to rear axle
reff: effective radius of the vehicle
wheel_speed: wheel speeds
1: FL
2: FR
3: RL
4: RR
'''
def get_slip_angles(self):
      slip_angle_1 = delta - (vy+lf*r)/vx*np.ones((2,))
      slip_angle_2 = (vy-lr*r)/vx*np.ones((2,))
      slip_angle = np.hstack((slip_angle_1,slip_angle_2))
      return slip_angle

def get_slip_ratio(self):
    slip_ratio = np.zeros((4,))
    for i in range(4):
        slip_ratio[i] = (reff*wheel_speed[i]-u)/(reff*wheel_speed[i])
    return slip_ratio