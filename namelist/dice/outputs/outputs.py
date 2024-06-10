#100 kpc Box with 13 levels of refinement from astropy import units as u
from astropy.constants import G, k_B, m_p
import astropy.units as u
from math import sqrt, pi
l_jeans = 12.22*u.pc*4
gamma = 5.0/3.0
Tgas = 100 * u.K
mu = 1.22
rho = pi * gamma *k_B * Tgas /(mu * G * l_jeans**2 * m_p)
rho_cgs = rho.decompose(u.cgs.bases)
rho_Msun_pc3 = rho_cgs.to(u.Msun/(u.pc**3))
print(rho.decompose(u.cgs.bases))
print(rho_Msun_pc3)
Mcell = rho_Msun_pc3 * (l_jeans/4)**3
print(Mcell)
nH = 0.75 * rho_cgs/m_p
print(nH.decompose(u.cgs.bases))

units_density = 0.677025430198932E-22 #1e9 Msol/kpc^3
units_time    = 0.470430312423675E+15 # G=1
units_length  = 0.308567758128200E+22 # kpc

units_mass = units_density * units_length**3
G_cgs = G.decompose(u.cgs.bases).value
print(G_cgs)
G_code = G_cgs * units_mass * units_time**2/units_length**3
print(G_code)
year_in_sec = u.year.to(u.second)
print(year_in_sec)
print("units_time_in Myr = ", units_time * 1.e-6/year_in_sec)
units_time_in_Myr = units_time * 1.e-6/year_in_sec
scale_length_in_cm = 1 * u.kpc.to(u.cm)
circle_len = scale_length_in_cm* 2 * pi
T_cgs = circle_len/(20*u.km.to(u.cm))
T_code = T_cgs/units_time
print(T_code)
print("period is ", T_code * units_time_in_Myr, " Myr")

print("delta_t_out = ", 10 * units_time_in_Myr)
print("1 Gyr in code units = ", 1000 * units_time_in_Myr)

