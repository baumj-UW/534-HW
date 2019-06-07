from numpy import pi, sqrt

# general inverter parameters
V = sqrt(2/3)*480
fg = 60
wg = 2*pi*fg
Tg = 1/fg
L = 100e-6 #100e-6
R = 0.75e-3 #0.75e-3
Vdc = 1250
tau = 1e-3
kp = L/tau
ki = R/tau

# per unit conversion base values
S_base = 1e6
w_base = wg
V_base = sqrt(2) * V
I = S_base / (3 * V) # Inom
I_base = sqrt(2) * I
Z_base = V / I

# per unit conversion
R = R / Z_base
ki = ki / (w_base * Z_base)
kp = kp / Z_base

# LCL filter
C0 = 100e-6 # TODO: find reasonable numbers for output cap
R0 = 0.1e-3
XL_base = w_base * L / (V / I)
B0 = w_base * C0 * (V / I)
R0 = R0 / Z_base