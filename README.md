# Roll-Yaw-and-Lateral-Velocity-MPC
Dallara AV21

This is for vehicle lateral stability of a general electirc vehicle as well as lateral position tracking of a Dallara AV-21 autonomous racing car.

Run the PYOMO modules on gooogle colab and the CasADi module on vs code after installing the required packages.

The vehicle dynamic model used is an in plane model with roll dynamics. For the Dallara AV-21 only  lateral control is considered. For the general EV we also consider the longitudinal control through an ABS system. Also note that the tire models used are linearized Pacejka and Dugoff tire models.


Some notation for reference:

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
