#----------------------------------------------------------------------------------#
# Import modules
import numpy as np

#----------------------------------------------------------------------------------#
# Define functions

def Q_AB(T_M, T_E, C_Pm, i_m, P_m, s, x):
    
    """The heat transferred from the melt to the cooling media 

    Keyword arguments:
    T_M -- Mold Temperature
    T_E -- Demolding Temperature
    C_Pm -- Specific heat of the melt 
    i_m -- Latent heat of fusion of polymer
    P_m -- Melt density
    s -- Part Thickness
    x -- Distance x 
    
	"""
    
    Q_AB = 10**-3 * ((T_M - T_E) * C_Pm + i_m) * (P_m * (s / 2) * x)
    return Q_AB

def a(k_m, P_m, C_Pm):
    
	"""Thermal diffusivity 

    Keyword arguments:
    k_m -- Thermal conductivity of the melt
    P_m -- Melt density
    C_Pm -- Specific heat of the melt 
    
    """
	a = k_m / (P_m * C_Pm)
	return a 

def Se(x, y, d):
    
	"""Shape Factor  

    Keyword arguments:
    x -- Distance x
    y -- Distance y
    d -- Diameter of cooling channel 
    
    """
	Se = (2 * np.pi) / (np.log(2 * x * np.sinh(((2 * np.pi * y) / x)) / (np.pi * d)))
	return Se

def Re(u, d, v):
    
	"""Reynolds number

    Keyword arguments:
    u -- Velocity of cooling water
    d -- Diameter of cooling channel 
    v -- Kinematic viscosity of cooling water 
    
    """
	Re = u * (d) / v
	return Re
 
def alpha(d, Re):
    
	"""Heat transfer coefficient 

    Keyword arguments:
    Re -- Reynolds number
    d -- Diameter of cooling channel 
    
    """
	alpha = (0.031395 / (d)) * (Re**0.8)
	return alpha

#----------------------------------------------------------------------------------#
# Basic Variables and Information

s = 2*10**-3           # Part thickness unit mm 
x = 30*10**-3        # Distance x unit mm
y = 10*10**-3         # Distance y unit mm
d = 10*10**-3          # Diameter of cooling channel unit mm
T_M = 250      # Mold temperature degree Celsius
T_E = 90      # Demolding temperature degree Celsius
i_m = 130      # Latent heat of fusion of polymer units kJ/kg
C_Pm = 2.5     # Specific heat of the melt units kJ/(kg*K)
P_m = 0.79     # Melt density units g/cm^3
k_m = 0.16     # Thermal conductivity of the melt units W/(m*K)
v = 1.2e-6     # Kinematic viscosity of cooling water units m^2/s
u = 1          # Velocity of cooling water units m/s
T_water = 15   # Temperature of cooling water degree Celsius
k_ST = 45      # Thermal conductivity of mold steel units W/(m*K)   
R= 5 
T_wall = (T_E + T_water) / 2 # Assumed mold wall temperature
print(T_wall)

#----------------------------------------------------------------------------------#
# Initial guess of wall temperature

iter = 1
err = 1.0
eps = 0.5 # Tolerance

while err > eps:
	# Calculation for the cooling time
	t_k = (s**2 / (np.pi**2 * a(k_m, P_m, C_Pm))) * \
		np.log((4 / np.pi) * ((T_M - T_wall) / (T_E - T_wall)))

	# Calculation for the heat received by the cooling agent (water)
	Q_w = (10**-3 * t_k) * (1 / (k_ST * Se(x, y, d)) + \
		(1 / (alpha(d, Re(u, d, v)) * (10**-3) * 2 * np.pi * R)))**-1 * (T_wall - T_water)
    
    # calculation for the mold temperature
	T_wall = T_water + (Q_w/((10**-3 * t_k) * (1 / (k_ST * Se(x, y, d)) + \
		(1 / (alpha(d, Re(u, d, v)) * (10**-3) * 2 * np.pi * R)))**-1)) # Guess 
	# Calculate err
	err = np.abs(Q_AB(T_M, T_E, C_Pm, i_m, P_m, s, x) - Q_w) / Q_AB(T_M, T_E, C_Pm, i_m, P_m, s, x)
	
	print('\nIteration: {}, err: {}'.format(iter, err))
	print("\n[INFO: Printing current values]")
	print('Q_w: {}'.format(Q_w))
	print('T_wall: {}'.format(T_wall))
	print('t_k: {}'.format(t_k))
	iter = iter + 1

#----------------------------------------------------------------------------------#
# Print data

print("\n\n[INFO: Iteration complete. Printing sample values]")
print("Q_AB: {} kJ/m".format(Q_AB(T_M, T_E, C_Pm, i_m, P_m, s, x)))
print ("a: {}".format(a(k_m, P_m, C_Pm)))
print ("Se: {}".format(Se(x, y, d)))
print("Reynold's: {}".format(Re(u, d, v)))
print("alpha: {} W/(m**2K)".format(alpha(d, Re(u, d, v))))

print("\n[INFO: Printing calculated values]")
print('T_wall:',T_wall) # Should be 37.83 C
print('t_k:',t_k) # Should be 8.03 s
print('Q_w:',Q_w) # Q_AB = Q_W

#----------------------------------------------------------------------------------#
# Plots
