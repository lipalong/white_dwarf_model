import math



#   typical core densities range from 10^8 to 10^12, in kg/m^3
#   typical radius of a solar mass white dwarf is on the order of 0.01 solar radii



#   number of electrons per nucleon; for an electrically neutral star, this is the same as the fraction of nucleons which are protons
#   for instance, for a star composed entirely of 12C (carbon) atoms, each nucleus has 6 protons and 6 neutrons

alpha = 0.5



#   mass of a nucleon, in kg

m_n = 1.67262178*10**(-27.0)



#   mass of an electron, in kg

m_e = 9.10938291*10**(-31.0)



#   mass of a hydrogen atom, in kg

m_H = 1.6605402*10**(-27.0)



#   Planck constant, in m^2*kg/s

h = 6.62606957*10**(-34.0)



#   gravitational constant, in m^3 kg^-1 s^-2

G = 6.67384*10**(-11.0)



#   speed of light, in m/s

c = 2.99792458*10**(8.0)



#   core density of the sun, in kg/m^3

rho_sun = 1.622*10**(5.0)



#   mass of the sun, in kg

m_sun = 1.98855*10**(30.0)



#   radius of the sun, in m

r_sun = 6.955*10**(8.0)



#   radius of the earth, in m

r_earth = 6.3781*10**(6.0)



#   K_F

K_F = ((3.0*0.5)/(8.0*math.pi*m_H))**(1.0/3.0)*(h/(m_e*c))



#   receive rho (core density) input

rho_input = input("Enter a value for the core density: ")



#   define scale factors

c_m = (4*math.pi*r_sun**(3.0)*rho_input)/m_sun
c_rho = (-1*G*m_H*m_sun)/(0.5*m_e*c**(2.0)*r_sun)
    
    
    
#   define Y(rho)

def Y(rho_input, rho, K_F):
    result = (K_F**(2.0)*rho_input**(2.0/3.0)*rho**(2.0/3.0))/(3*(1+K_F**(2.0)*rho_input**(2.0/3.0)*rho**(2.0/3.0))**(1.0/2.0))
    return result



#   define stellar structure equations

def mf(c_m, r, rho):
    result = c_m*r**(2.0)*rho
    return result

def rhof(c_rho, r, m, rho, Y):
    result = (c_rho*m*rho)/(r**(2.0)*Y(rho_input, rho, K_F))
    return result



#   begin euler's algorithm

h = 0.000000001
rho = rho_input/rho_sun
r = h
m = 4.0*math.pi*r**2*rho

print "{}\t{}\t{}".format("radius", "mass", "density")

while rho>=0:
    print "{}\t{}\t{}".format(r, m, rho)
    r += h
    m += h*mf(c_m, r, rho)
    rho += h*rhof(c_rho, r, m, rho, Y)
    
print "{}\t{}\t{}".format(r, m, rho)



#   scaling the radius to earth radii

r_es = r*(r_sun/r_earth)



#   print results

print "for a core density of",rho_input,"kg/m^3, the radius is",r_es, "earth radii and the mass is",m, "solar masses."
