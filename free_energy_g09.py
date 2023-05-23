"""

 Thermochemistry based on gaussian09 FREQ calculation

 python3 free_energy_g09.py [g09file.log] [temperature]

"""

import numpy as np
import sys
import readLog

R = 0.00000316679085237                                     # Kelvin to Hartree
NA = 6.02214086e23
kB = R/NA
P_conv = 3.4439667564003e-9                                 # atm to a.u.
mass_conv = 1836.15267389                                   # a.m.u to a.u.
vib_temp_conv = 1.438786296                                 # cm-1 to (equivalent) Kelvin
anharmonic_corr = 0.9854                                    # empirical anharmonic correction for the frequencies.
#anharmonic_corr = 1.0                                    # empirical anharmonic correction for the frequencies.


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
    print(" (T = " + str(T) + ")") 
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

    LogFile = sys.argv[1]
    hp = readLog.catch_hp_g09(LogFile)
    freq = readLog.catch_freq_g09(LogFile, hp)
    vib_temp = anharmonic_corr*vib_temp_conv*np.array(freq) # vibrational temperatures

    T = float(sys.argv[2])                                  # Temperature (Kelvin)
    print(" => Vibrational analysis at T = %f K <= \n"%T)
    P = 1.0                                                 # Pressure (atm)
    P *= P_conv

    spin_mult = readLog.catch_spin_mult(LogFile)            # spin multiplicity
    Ix, Iy, Iz = readLog.catch_rot_temp(LogFile)            # Rotational temperatures (Kelvin)
    sigmar = readLog.catch_rot_symmetry_number(LogFile)     # Rotational symmetry number
    mass = readLog.catch_mol_mass(LogFile)                  # mass (a.m.u)
    mass *= mass_conv

    E0 = Etot(0., vib_temp)
    E = Etot(T, vib_temp)
    S = Stot(mass, T, P, Ix, Iy, Iz, sigmar, vib_temp, spin_mult)
    print(" Free energy correction at T = 0K (ZPE):\t" + str(E0))
    print(" Internal energy correction at T = %fK:\t"%T + str(E))
    print(" Free energy correction at T = %fK:\t"%T + str(E+T*R-T*S))
