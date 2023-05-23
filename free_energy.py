"""

 Thermochemistry based on gaussian09 FREQ calculation

 python3 free_energy.py [freq_file.dat] [temperature]

"""

import numpy as np
import sys

R = 0.00000316679085237                                     # Kelvin to Hartree
NA = 6.02214086e23
kB = R/NA
P_conv = 3.4439667564003e-9                                 # atm to a.u.
mass_conv = 1836.15267389                                   # a.m.u to a.u.

# Energy
def Et(T):
    return 1.5*R*T

def Er(T):
    return 1.5*R*T

def Ev(T, vib_temp):
    nvib = len(vib_temp)
    Ev = 0 
    if T == 0.:
        for i in range(nvib):
            theta = vib_temp[i]
            Ev += theta*0.5
        return R*Ev
    for i in range(nvib):
        theta = vib_temp[i]
        Ev += theta*(0.5 + 1/(np.exp(theta/T)-1))
    return R*Ev
 
def Etot(T, vib_temp):
    print(" Translational energy:\t" + str(Et(T)))
    print(" Rotational energy:\t" + str(Er(T)))
    print(" Vibrational energy:\t" + str(Ev(T,vib_temp)) + "\n")

    return Et(T) + Er(T) + Ev(T,vib_temp)

# Entropy
def St(mass, T, P):
    Zt = (mass*R*T/(2*np.pi))**1.5
    qt = Zt*R*T/P
    return R*(np.log(qt) + 2.5)

def Sr(Ix, Iy, Iz, T, sigmar):
    Zr = T**1.5/np.sqrt(Ix*Iy*Iz)
    qr = np.sqrt(np.pi)*Zr/sigmar
    return R*(np.log(qr) + 1.5)

def Sv(T, vib_temp):
    nvib = len(vib_temp)
    Sv = 0
    for i in range(nvib):
        theta = vib_temp[i]/T
        Sv += theta/(np.exp(theta)-1) - np.log(1-np.exp(-theta))
    return R*Sv

def Se(spin_mult):
    return R*np.log(spin_mult)

def Stot(mass, T, P, Ix, Iy, Iz, sigmar, vib_temp, spin_mult):
    print(" Translational entropy:\t" + str(St(mass,T,P)))
    print(" Rotational entropy:\t" + str(Sr(Ix,Iy,Iz,T,sigmar)))
    print(" Vibrational entropy:\t" + str(Sv(T,vib_temp)))
    print(" Electronic entropy:\t" + str(Se(spin_mult)) + "\n")

    return St(mass, T, P) + Sr(Ix, Iy, Iz, T, sigmar) + Sv(T, vib_temp) + Se(spin_mult)

if __name__=="__main__":

    f = open(sys.argv[1])
    vib_temp = [0.9854*float(line) for line in f.readlines()]      # vibrational temperatures
    f.close()

    T = float(sys.argv[2])                                  # Temperature (Kelvin)
    print(" => Vibrational analysis at T = %f K <= \n"%T)
    P = 1.0                                                 # Pressure (atm)
    P *= P_conv

    mass = 210.08521                                        # mass (a.m.u)
    spin_mult = 1                                           # spin multiplicity
    Ix, Iy, Iz = 0.04329, 0.02375, 0.01743                  # Rotational temperatures (Kelvin)
    sigmar = 1                                              # Rotational symmetry number

    mass *= mass_conv

    E0 = Etot(0., vib_temp)
    E = Etot(T, vib_temp)
    S = Stot(mass, T, P, Ix, Iy, Iz, sigmar, vib_temp, spin_mult)
    print(" Vib. free energy at T = 0K (ZPE):\t" + str(E0))
    print(" Vib. energy at T = %fK:\t"%T + str(E))
    print(" Vib. free energy at T = %fK:\t"%T + str(E+T*R-T*S))