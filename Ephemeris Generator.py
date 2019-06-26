# Ephemeris Generator

from visual import *
from math import *
from numpy import *

def Eiteration(e,M,Ei): #to calculate E by iteration
    E=Ei-((M-(Ei-e*sin(Ei)))/(e*cos(Ei)-1))
    if abs(M-(E-e*sin(E))<=1*10**(-9)):
        return E
    else:
        return Eiteration(e,M,E)

def ephemerisGen(orbElm, JDfloat, SunVectorArray):
    #defines the elements in orbElm and assigns them positions in the array
    a=orbElm[0]; e=orbElm[1]; i=orbElm[2]*pi/180; omega=orbElm[3]*pi/180; w=orbElm[4]*pi/180; T=orbElm[5]
    #gives them values
##    a=0.9223140
##    e=0.1910619
##    i=.0581541961 #rad-i can mod this 
##    omega=3.567933131 #rad
##    w=2.20642482 #rad
##    T=2456189.09508 #JD
    #defines the Sun vecor components:
    Rx=SunVectorArray[0]; Ry=SunVectorArray[1]; Rz=SunVectorArray[2]
    #gives them values:
##    Rx=0.91740051 
##    Ry=0.371820921
##    Rz=0.161116271
    R=matrix([[Rx],[Ry],[Rz]])
##    JDfloat=2462240.416666667
    
    k=.01720209894
    n=float(k*((1.000/(a**3))**.5))
    M=float(n*(JDfloat-T))

    E0=M
    E=Eiteration(e,M,E0)

    x=a*cos(E)-a*e
    y=a*(1-e**2)**.5*sin(E)
    z=0
    cartesian=matrix([[x],[y],[z]])
    #cartesian coords of pos vector
##    print cartesian

    #first transform
    warray=matrix([[cos(w), (-1)*sin(w), 0],[sin(w), cos(w), 0],[0,0,1]])
    mult_1=warray*cartesian
##    return mult_1
    
    #second transform
    iarray=matrix([[1,0,0],[0,cos(i),(-1)*sin(1)],[0,sin(i),cos(i)]])
    mult_2=iarray*mult_1
##    return mult_2

    #ecliptic x1, y1, z1
    omegaarray=matrix([[cos(omega),(-1)*sin(omega),0],[sin(omega),cos(omega),0],[0,0,1]])
    ecliptic=omegaarray*mult_2
##    print ecliptic

    fancyE=23.4376600557*pi/180 #obliquity
    #equatorial x2, y2, z2
    Earray=matrix([[1,0,0],[0,cos(fancyE),(-1)*sin(fancyE)],[0,sin(fancyE),cos(fancyE)]])
    r=Earray*ecliptic
    rmag=linalg.norm(r)
##    print r

    #now to get the rho vector
    rho=r+R
##    print rho
    rhoDist=linalg.norm(rho) #magnitude
    rhoHat=rho/rhoDist #unit vector
    rhoHatx=rhoHat[0]
    rhoHaty=rhoHat[1]
    rhoHatz=rhoHat[2]

    #to calculate RA and dec
    dec=asin(rhoHatz)
    RAsin1 = asin(rhoHaty/cos(dec))
    RAcos1 = acos(rhoHatx/cos(dec))
    RAcos2 = -RAcos1
##    print RAsin1, RAcos1, RAcos2
    if round(RAsin1,11) == round(RAcos2,11):
        RA = min(RAsin1,RAcos2)
    else:
        RA = abs(RAsin1)
    decDeg=dec*180/pi
##    print RA
    RAhours=((RA)*(180/pi)/15) %24
##    print RAhours, decDeg

    solution=[RAhours,decDeg,rhoDist]
    return solution

#################################################

orbElm = [2.06087341461,0.823506396318,0.276269946,4.496898739,4.329961048,0.744206218]
JDfloat = 2456107.897319
SunVectorArray = [-0.1392320229883591,0.9239200060886814,0.4005503954711760]

print ephemerisGen(orbElm,JDfloat,SunVectorArray)

# solutions to a:[20.18324820681805, -18.578553721735656]--[RA,dec]
# solution to b: 0.020457246246584442 AU-this is when it's supposed to be crashing into the earth
# solutions to c: [20.6, -15.6, 0.0571]
# they could have downgraded the risk using more sig figs
