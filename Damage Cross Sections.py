# -*- coding: utf-8 -*-
# This code was written by Barnaby Levin of the Muller Group at Cornell University.
# It uses formulae from a report by CR Bradley of Argonne National Labs to compute displacement damage cross sections
# for atoms subjected to electron radiation. 
# The data will be saved in text format, and a preview plot will be displayed. 

import tkinter
import math
import numpy as np
import matplotlib.pyplot as plt

Name = 0

Window = tkinter.Tk()

#Set title, size, and colour of window
Window.wm_title("User Input")
Window.geometry("400x620")
Window.resizable(1000,1000)
Window.configure(background='black')

# Code to add widgets goes here
# Below is code to allow the user to enter text. 
# The code for each text box is labelled with a number.  
# Data to input can be found in Bradley's report, or from other sources

#Sp = Tkinter.Label(Window,text=" ", bg="black",)
#Sp.pack()
L0 = tkinter.Label(Window,text="Please Enter Info Below", font="Helvetica 20 bold", bd=10, bg="black", fg = "white")
L0.pack()
Sp0 = tkinter.Label(Window,text=" ", bg="black",)
Sp0.pack()

#This is the name that files will be saved with
L1 = tkinter.Label(Window,text="Substance Name: ", font="Helvetica 12 bold", bd=5, bg="black", fg = "white")
L1.pack()
E1 = tkinter.Entry(Window) 
E1.pack()
E1.insert(0, "Enter Substance Name Here")
Sp1 = tkinter.Label(Window,text=" ", bg="black",)
Sp1.pack()

#The next few inputs are for the calculations
L2 = tkinter.Label(Window,text="Atomic Number Z: ", font="Helvetica 12 bold", bd=5, bg="black", fg = "white")
L2.pack()
E2 = tkinter.Entry(Window) 
E2.pack()
E2.insert(0, "3")
Sp2 = tkinter.Label(Window,text=" ", bg="black")
Sp2.pack()

L3 = tkinter.Label(Window,text="Atomic Mass in amu: ", font="Helvetica 12 bold", bd=5, bg="black", fg = "white")
L3.pack()
E3 = tkinter.Entry(Window) 
E3.pack()
E3.insert(0, "6.941")
Sp3 = tkinter.Label(Window,text=" ", bg="black")
Sp3.pack()

#L4 = Tkinter.Label(Window,text="Density of material (g/cm^3): ",bd=5, bg="black", fg = "white")
#L4.pack()
#E4 = Tkinter.Entry(Window) 
#E4.pack()
#E4.insert(0, "0.534")
#Sp4 = Tkinter.Label(Window,text=" ", bg="black")
#Sp4.pack()

L5 = tkinter.Label(Window,text="Displacement Energy (eV): ", font="Helvetica 12 bold", bd=5, bg="black", fg = "white")
L5.pack()
E5 = tkinter.Entry(Window) 
E5.pack()
E5.insert(0, "2.54")
Sp5 = tkinter.Label(Window,text=" ", bg="black")
Sp5.pack()

L6 = tkinter.Label(Window,text="Voltage Step Size For Calculations (keV): ", font="Helvetica 12 bold", bd=5, bg="black", fg = "white")
L6.pack()
E6 = tkinter.Entry(Window) 
E6.pack()
E6.insert(0, "0.5")
Sp6 = tkinter.Label(Window,text=" ", bg="black")
Sp6.pack()
Sp7 = tkinter.Label(Window,text=" ", bg="black")
Sp7.pack()

L7 = tkinter.Label(Window,text="Damage Type: ", font="Helvetica 12 bold", bd=5, bg="black", fg = "white")
L7.pack()

# The code will execute when the user clicks the button for the type of damage they want to calculate for.
# Let us define the function that will be executed for VED and surface sputtering.
# This will use the McKinley Feshbach formula to calculate damage cross section.
# First task - define the function that happens on the button click i.e. the calculations
def callback1():
    
    #First part of the function grabs data from the window, coverts it to the data type we want
    Name=E1.get()
    Ans2=E2.get()
    Ans3=E3.get()
#    Ans4=E4.get()
    Ans5=E5.get()
    Ans6=E6.get()
    print("Calculating......")
    Savename = str(Name) + ".txt"
    Z=float(Ans2)
    M=float(Ans3)
 #   Dens=float(Ans4)
    Ed=float(Ans5)
    No=float(Ans6)

    
    #Now define a few constants for our calculations
    q = 1.6*(10**-19)
    M0 = M*931494000 # Mass of nucleus in eV
    m0 = 511000 # Mass of electron in eV
    pi = 3.14159265359
    hbar = 1.05457148*(10**-34)
#    h = 4.1356675*(10**-15)
    c = 299792458
    epsilon0 = 8.854187818*(10**-12)
 #   k = 1/(4*pi*epsilon0)
  #  Na = (10**6)*Dens/(M*1.66*(10**-24)) # Number density of atoms
    alpha = Z*q*q/(4*pi*epsilon0*hbar*c)
#    a0 = 5.2917721*(10**-11) 
#    R = a0/(Z**(1/3))
    N = int(500/No)
    
    #Now define some helpful arrays for storing data
    EkeV = np.ones(N)
    EeV = np.ones(N)
    SigmaMF = np.ones(N)
    
    
    #OK, let's fill up the arrays - for loop time!
    for n1 in range (0,N):
        EkeV[n1] = n1*No
        EeV[n1] = 1000*n1*No
        
        # Great, now let's start computing cross sections
        v = math.sqrt(1-((1+(EeV[n1]/m0))**-2))
        Tm = 2*EeV[n1]*(EeV[n1]+2*m0)/M0
        
        #McKinley - Fesbach using formulae adapted from Bradley report for Argonne NL
        if Ed > Tm:
                SigmaMF[n1] = 0;
        else:
                f0 = (Z*q/(m0*4*pi*epsilon0))**2
                f1 = 2*(math.sqrt(Tm/Ed)-1)-math.log(Tm/Ed)
                Zeta = (((Tm/Ed)-1)-((v**2)*math.log(Tm/Ed))+pi*alpha*v*f1)
                f2 = pi*f0*(1-(v**2))/(v**4)
                SigmaMF[n1] = f2*Zeta*(10**28)
                if(SigmaMF[n1] < 0):
                    SigmaMF[n1] = 0
                                       
    #Now that we have the cross sections, python will output a preview figure
    
    plt.figure(figsize=(12,8), facecolor='w') #Defines figure with white background, 12'' wide by 8'' tall
    Plot = plt.plot(EkeV, SigmaMF) #Plots MF and KP cross sections as function of beam voltage
    plt.rc('axes', linewidth=2) # Sets width of axes
    plt.rcParams['xtick.major.size'] = 6 #Sets height of x axis ticks
    plt.rcParams['xtick.major.width'] = 2 #Sets width of x axis ticks
    plt.rcParams['ytick.major.size'] = 6 #Same as above but for y axis ticks
    plt.rcParams['ytick.major.width'] = 2
    plt.xticks(size=20, family='Arial') #Set font size and font type for ticks
    plt.yticks(size=20, family='Arial')
    
    # Set the title and axis labels with sizes, fonts and font weights
    plt.title('\nKnock-on displacement damage cross section\nas a function of electron beam voltage',size=20, fontweight='bold', family='Arial')
    plt.xlabel('Electron Beam Voltage (kV)', size=20, fontweight='bold', family='Helvetica')
    plt.ylabel('Cross Section (barns)', size=20, fontweight='bold', family='Helvetica')
    
    plt.setp(Plot, color='r', linewidth=2.0)
        # Make the plot appear
    plt.show()
    
    # Now create a storage array for writing data to txt file
    Store = np.ones((N,2))
    for n2 in range (0,N): #Put all the data in the file
        Store[n2,0] = EkeV[n2]
        Store[n2,1] = SigmaMF[n2]
    
    # Write the data                                                                                                                                       
    np.savetxt(Savename,Store,delimiter='\t',newline='\n')
    print(Savename, "..........saving")
    Window.destroy() #Closes window

# Now we define the command used for bulk displacement calculations. 
# This will use the Kinchin Pease formula for damage to account for the possibility of atomic cascades. 
# First task - define the function that happens on the button click i.e. the calculations
def callback2():
    
    #First part of the function grabs data from the window, coverts it to the data type we want
    Name=E1.get()
    Ans2=E2.get()
    Ans3=E3.get()
#    Ans4=E4.get()
    Ans5=E5.get()
    Ans6=E6.get()
    print("Calculating......")
    Savename = str(Name) + ".txt"
    Z=float(Ans2)
    M=float(Ans3)
 #   Dens=float(Ans4)
    Ed=float(Ans5)
    No=float(Ans6)

    
    #Now define a few constants for our calculations
    q = 1.6*(10**-19)
    M0 = M*931494000 # Mass of nucleus in eV
    m0 = 511000 # Mass of electron in eV
    pi = 3.14159265359
    hbar = 1.05457148*(10**-34)
#    h = 4.1356675*(10**-15)
    c = 299792458
    epsilon0 = 8.854187818*(10**-12)
 #   k = 1/(4*pi*epsilon0)
  #  Na = (10**6)*Dens/(M*1.66*(10**-24)) # Number density of atoms
    alpha = Z*q*q/(4*pi*epsilon0*hbar*c)
#    a0 = 5.2917721*(10**-11) 
#    R = a0/(Z**(1/3))
    N = int(500/No)
    
    #Now define some helpful arrays for storing data
    EkeV = np.ones(N)
    EeV = np.ones(N)
    SigmaMF = np.ones(N)
    SigmaKP = np.ones(N)
    
    
    #OK, let's fill up the arrays - for loop time!
    for n1 in range (0,N):
        EkeV[n1] = n1*No
        EeV[n1] = 1000*n1*No
        
        # Great, now let's start computing cross sections
        v = math.sqrt(1-((1+(EeV[n1]/m0))**-2))
        Tm = 2*EeV[n1]*(EeV[n1]+2*m0)/M0
        
        #McKinley - Fesbach first, Kinchin Pease must be greater than or equal to this. 
        if Ed > Tm:
                SigmaMF[n1] = 0;
        else:
                f0 = (Z*q/(m0*4*pi*epsilon0))**2
                f1 = 2*(math.sqrt(Tm/Ed)-1)-math.log(Tm/Ed)
                Zeta = (((Tm/Ed)-1)-((v**2)*math.log(Tm/Ed))+pi*alpha*v*f1)
                f2 = pi*f0*(1-(v**2))/(v**4)
                SigmaMF[n1] = f2*Zeta*(10**28)
                if(SigmaMF[n1] < 0):
                    SigmaMF[n1] = 0
                        
        # Now Kinchin - Pease, the cross section accounting for cascades.
        if(Ed>Tm):
                SigmaKP[n1] = 0
        else:
                F0 = pi*((Z*q/(m0*4*pi*epsilon0))**2)*((1-v*v)/(v*v*v*v))
                F1 = (Tm/(2*Ed)) - (pi*alpha*v+v*v)*math.log(2) + 2*pi*alpha*v*(math.sqrt(Tm/Ed) - math.sqrt(Tm/(2*Ed)))
                F2 = math.log(Tm/(2*Ed)) + v*v*(((2*Ed)/Tm)-1) + pi*alpha*v*(((2*Ed)/Tm) - 2*math.sqrt((2*Ed)/Tm) + 1)
                SigmaKP[n1] = F0*(F1 + 0.5*(Tm/Ed)*F2)*(10**28)
                if(SigmaKP[n1]<SigmaMF[n1]):
                    SigmaKP[n1]=SigmaMF[n1]                
    
    #Now that we have the cross sections, python will output a preview figure
    
    plt.figure(figsize=(12,8), facecolor='w') #Defines figure with white background, 12'' wide by 8'' tall
    Plot = plt.plot(EkeV, SigmaKP) #Plots MF and KP cross sections as function of beam voltage
    plt.rc('axes', linewidth=2) # Sets width of axes
    plt.rcParams['xtick.major.size'] = 6 #Sets height of x axis ticks
    plt.rcParams['xtick.major.width'] = 2 #Sets width of x axis ticks
    plt.rcParams['ytick.major.size'] = 6 #Same as above but for y axis ticks
    plt.rcParams['ytick.major.width'] = 2
    plt.xticks(size=20, family='Arial') #Set font size and font type for ticks
    plt.yticks(size=20, family='Arial')
    
    # Set the title and axis labels with sizes, fonts and font weights
    plt.title('\nKnock-on displacement damage cross section\nas a function of electron beam voltage',size=20, fontweight='bold', family='Arial')
    plt.xlabel('Electron Beam Voltage (kV)', size=20, fontweight='bold', family='Arial')
    plt.ylabel('Cross Section (barns)', size=20, fontweight='bold', family='Arial')

    # Define line colour and line width on the plot
    plt.setp(Plot, color='r', linewidth=2.0)
    # Make the plot appear
    plt.show()
    
    # Now create a storage array for writing data to txt file
    Store = np.ones((N,2))
    for n2 in range (0,N): #Put all the data in the file
       Store[n2,0] = EkeV[n2]
       Store[n2,1] = SigmaKP[n2]   
    
    # Write the data                                                                                                                                       
    np.savetxt(Savename,Store,delimiter='\t',newline='\n')
    print(Savename, "..........saving")
    Window.destroy() #Closes window

# Now define the button. Set to execute the command defined above
B1 = tkinter.Button(Window, text="Vacancy Enhanced Displacement", width=30, command=callback1)
B1.pack() 
Sp8 = tkinter.Label(Window,text=" ", bg="black")
Sp8.pack()
B2 = tkinter.Button(Window, text="Surface Sputtering", width=20, command=callback1)
B2.pack() 
Sp9 = tkinter.Label(Window,text=" ", bg="black")
Sp9.pack()
B3 = tkinter.Button(Window, text="Bulk Displacement", width=20, command=callback2)
B3.pack() 
Sp9 = tkinter.Label(Window,text=" ", bg="black")
Sp9.pack()

# Tell Python this is the end of the window code 
Window.mainloop()

