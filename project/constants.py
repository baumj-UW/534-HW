from numpy import pi, sqrt

# general inverter parameters
V = sqrt(2/3)*480
fg = 60
wg = 2*pi*fg
Tg = 1/fg
Vdc = 1250

# LCL filter
L = 10e-6 #100e-6
R = 0.007
Lg = 2e-6
Rg = 0.0012
C0 = 2.4e-3 # TODO: find reasonable numbers for output cap
R0 = 1

# current controller
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

# per unit controller
ki = ki / (w_base * Z_base)
kp = kp / Z_base

# per unit LCL
XL = w_base * L / (V / I)
XLg = w_base * Lg / (V / I)
R = R / Z_base
Rg = Rg / Z_base
B0 = w_base * C0 * (V / I)
R0 = R0 / Z_base