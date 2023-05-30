syms ca slip_angle slip_ratio mu Fz cs
lamb = mu*Fz*(1+slip_ratio)/(2*((cs*slip_ratio)^2+(ca*tan(slip_angle)^2)^0.5));
flamb = lamb*(2-lamb)
Fy = ca*tan(slip_angle)*flamb/(1+slip_ratio)
diff(Fy,slip_angle)