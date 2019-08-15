# -*- coding: utf-8 -*-
# This code was written by Barnaby D.A. Levin of the Muller Group at Cornell University.
# The code uses the McKinley-Feshbach Formula to calculate the number of knock-on damage events per atom for a given set of imaging parameters in an electron microscope.  
# The code requires activation energy for atomic displacement to be known. These can be found in literature. 
"""
Created on Fri Jul 10 13:37:28 2015

@author: Barnaby Levin
"""
import math

Z=float(input('Enter Atomic Number Z:\t\n'))
M=float(input('Enter Atomic Mass in A.M.U.:\t\n'))
kV=float(input('Enter Beam Voltage in kV:\t\n'))
I=float(input('Enter Beam Currernt in pA:\t\n'))
FOV=float(input('Enter Image Field of View in nm:\t\n'))
Pix=float(input('Enter Image Width in Pixels:\t\n'))
DwT=float(input('Enter Dwell Time per Pixel in us:\t\n'))
No=float(input('Enter Number of Images Acquired:\t\n'))
Ed=float(input('Enter Atomic Displacement Energy in eV:\t\n')) 

#Now define a few constants for our calculations
q = 1.6*(10**-19)
V=kV*1000
M0 = M*931494000 # Mass of nucleus in eV
m0 = 511000 # Mass of electron in eV
pi = 3.14159265359
hbar = 1.05457148*(10**-34)
#    h = 4.1356675*(10**-15)
c = 299792458
epsilon0 = 8.854187818*(10**-12)
alpha = Z*q*q/(4*pi*epsilon0*hbar*c)
v = math.sqrt(1-((1+(V/m0))**-2))
Tm = 2*V*(V+2*m0)/M0
if Ed > Tm:
        SigmaMF = 0;
else:
        f0 = (Z*q/(m0*4*pi*epsilon0))**2
        f1 = 2*(math.sqrt(Tm/Ed)-1)-math.log(Tm/Ed)
        Zeta = (((Tm/Ed)-1)-((v**2)*math.log(Tm/Ed))+pi*alpha*v*f1)
        f2 = pi*f0*(1-(v**2))/(v**4)
        SigmaMF = f2*Zeta*(10**28)
        if(SigmaMF < 0):
            SigmaMF = 0
            
#Now calculate Dose per nm2.
Dose=I*0.000000000001*DwT*Pix*Pix*0.000001*No/(FOV*FOV*q)
Damage=Dose*SigmaMF*0.0000000001
print('Knock-on Damage Cross Section (barns): ', SigmaMF)
print('Total Dose (electrons per nm2): ', Dose)
print('Damage events per target atom: ', Damage)
