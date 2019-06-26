from __future__ import division
from math import *
from visual import *
import numpy
import csv

####################################################################################################

# Orbit Determination, 20 July 2012
# Written by John Khouri
# The purpose of this program is to generate the six orbital elements from the RA and declination
# of an asteroid at three given times.

####################################################################################################
####################################################################################################

def degs2rads(degs):
    # Convert degrees to radians
    rads = (degs/360)*2*pi
    return rads
    
def hrs2rads(hrs):
    # Convert hours to radians
    rads = (hrs*pi)/12
    return rads
    
def dot(a,b):
    # Take the dot product of two vectors
    solution = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    return solution

def cross(a,b):
    # Take the cross product of two vectors
    solution = [a[1]*b[2]-a[2]*b[1],-a[0]*b[2]+a[2]*b[0],a[0]*b[1]-a[1]*b[0]]
    return solution

def JulianDate(month,day,year,time):
    J0 = 367*year - int((7*(year + int((month+9)/12)))/4) + int((275*month)/9) + day + 1721013.5
    JD = J0 + time/24
    return JD

def UnitVector(RA,Dec):
    # Return the unit vector of a point given its RA and Dec
    p_Hat = [cos(RA)*cos(Dec),sin(RA)*cos(Dec),sin(Dec)]
    return p_Hat

def Vector(x,RA,Dec):
    # Use magnitude (x), RA, and Dec to find a vector
    y = UnitVector(RA,Dec)
    xvec = [x*y[0],x*y[1],x*y[2]]
    return xvec

def f(rd,r2,rdot,tau):
    f = 1 - (tau**2)/(2*rd**3) + ((tau**3)*dot(r2,rdot))/(2*(rd**5)) + ((tau**4)/24)*((3/(rd**3))*((dot(rdot,rdot)/(rd**2))-1/(rd**3))-(15/(rd**2))*dot(rdot,rdot)**2+1/(rd**6))
    return f

def g(rd,rdot,tau):
    g = tau - (tau**3)/(6*rd**3) + ((tau**4)*dot(rdot,rdot))/(4*(rd**5))
    return g

def a(rd,r2,rdot,tau1,tau3):
    # Reapproximate a1 and a3
    a1 = (g(rd,rdot,tau3))/(f(rd,r2,rdot,tau1)*g(rd,rdot,tau3) - f(rd,r2,rdot,tau3)*g(rd,rdot,tau1))
    a3 = -(g(rd,rdot,tau1))/(f(rd,r2,rdot,tau1)*g(rd,rdot,tau3) - f(rd,r2,rdot,tau3)*g(rd,rdot,tau1))
    return a1,a3

def VectorRDot(r1,r2,r3,rd,rdot,tau1,tau3):
    # Calculate r dot
    x = f(rd,r2,rdot,tau3)
    print x
    y = g(rd,rdot,tau3)
    print y
    r1dot = [(r3[0]-(r2[0]*x))/y,(r3[1]-(r2[1]*x))/y,(r3[2]-(r2[2]*x))/y]
    print rdot
    x = f(rd,r2,rdot,tau1)
    y = g(rd,rdot,tau1)
    r2dot = [(r1[0]-(r2[0]*x))/y,(r1[1]-(r2[1]*x))/y,(r1[2]-(r2[2]*x))/y]
    rdot = [(r1dot[0] + r2dot[0])/2,(r1dot[1] + r2dot[1])/2,(r1dot[2] + r2dot[2])/2]
    return rdot

def AngleAmb(angle1,angle2):
    if angle1 < 0:
            angle1m = [pi - angle1, 2*pi + angle1]
    else:
            angle1m = [angle1, pi - angle1]
    if angle2 < 0:
        angle2m = [-angle2, 2*pi + angle2]
    else:
        angle2m = [angle2, 2*pi - angle2]
    for angle in angle2m:
        if abs(angle-angle1m[0]) < 1*10**-12:
            return angle
        elif abs(angle-angle1m[1])<1*10**-12:
            return angle
    return angle1m,angle2m
    
####################################################################################################

def Recursion(RA1,Dec1,t1,RA2,Dec2,t2,RA3,Dec3,t3,R1,R2,R3,R1dot,R2dot,R3dot):

    # Convert RA and Dec into radians (Correct)
    RA1 = hrs2rads(RA1)
    RA2 = hrs2rads(RA2)
    RA3 = hrs2rads(RA3)
    Dec1 = degs2rads(Dec1)
    Dec2 = degs2rads(Dec2)
    Dec3 = degs2rads(Dec3)
    
    # Define k
    k = 0.01720209895

    
    
    # Calculate tau (Correct)
    tau1 = k*(t1-t2)
    tau2 = k*(t3-t1)
    tau3 = k*(t3-t2)

    # Calculate unit vector p
    #p1_Hat = UnitVector(RA1,Dec1) (Correct)
    p1_Hat = [cos(RA1)*cos(Dec1),sin(RA1)*cos(Dec1),sin(Dec1)]
    p2_Hat = [cos(RA2)*cos(Dec2),sin(RA2)*cos(Dec2),sin(Dec2)]
    p3_Hat = [cos(RA3)*cos(Dec3),sin(RA3)*cos(Dec3),sin(Dec3)]
    
    # Estimate a1 and a3 (Correct)
    a1 = tau3/tau2
    a3 = -tau1/tau2

    # Calculate the magnitude of each p, pd
    p1d = (a1*dot(cross(R1,p2_Hat),p3_Hat) - dot(cross(R2,p2_Hat),p3_Hat) + a3*dot(cross(R3,p2_Hat),p3_Hat))/(a1*dot(cross(p1_Hat,p2_Hat),p3_Hat))
    p2d = (a1*dot(cross(p1_Hat,R1),p3_Hat) - dot(cross(p1_Hat,R2),p3_Hat) + a3*dot(cross(p1_Hat,R3),p3_Hat))/(-1*dot(cross(p1_Hat,p2_Hat),p3_Hat))
    p3d = (a1*dot(cross(p2_Hat,R1),p1_Hat) - dot(cross(p2_Hat,R2),p1_Hat) + a3*dot(cross(p2_Hat,R3),p1_Hat))/(a3*dot(cross(p2_Hat,p3_Hat),p1_Hat))

    # Correct for light speed change
    c = 173.1446 # speed of light in AU per day(JD)
    t1 = t1 - p1d/c
    t2 = t2 - p2d/c
    t3 = t3 - p3d/c

    # Recalculate tau after light speed change
    tau1 = k*(t1-t2)
    tau2 = k*(t3-t1)
    tau3 = k*(t3-t2)
    print "tau =", tau1,tau2,tau3

    
    # Calculate vector p 
    p1 = Vector(p1d,RA1,Dec1)
    p2 = Vector(p2d,RA2,Dec2)
    p3 = Vector(p3d,RA3,Dec3)
    
    # Calculate vector r and magnitude of r2
    r1 = [p1[0]-R1[0],p1[1]-R1[1],p1[2]-R1[2]]
    r2 = [p2[0]-R2[0],p2[1]-R2[1],p2[2]-R2[2]]
    r3 = [p3[0]-R3[0],p3[1]-R3[1],p3[2]-R3[2]]
    rd = mag(r2)
         
    # Use linear approximation to approximate rdot
    rdot = [(r3[0]-r1[0])/(tau3-tau1),(r3[1]-r1[1])/(tau3-tau1),(r3[2]-r1[2])/(tau3-tau1)]

    # Calculate f1,g1,f3, and g3
    f1 = f(rd,r2,rdot,tau1)
    g1 = g(rd,rdot,tau1)
    f3 = f(rd,r2,rdot,tau3)
    g3 = g(rd,rdot,tau3)

    print "f1 =",f1
    print "g1 =",g1
    print "f3 =",f3
    print "g3 =",g3

    

    print "r =", r2
    print "rdot =", rdot


    ## Second iteration
    print " "
    print " "

    
    # Recalculate a1 and a3
    a1,a3 = a(rd,r2,rdot,tau1,tau3)
    
    print "a1,a3 =", a1,a3

    p1d = (a1*dot(cross(R1,p2_Hat),p3_Hat) - dot(cross(R2,p2_Hat),p3_Hat) + a3*dot(cross(R3,p2_Hat),p3_Hat))/(a1*dot(cross(p1_Hat,p2_Hat),p3_Hat))
    p2d = (a1*dot(cross(p1_Hat,R1),p3_Hat) - dot(cross(p1_Hat,R2),p3_Hat) + a3*dot(cross(p1_Hat,R3),p3_Hat))/(-1*dot(cross(p1_Hat,p2_Hat),p3_Hat))
    p3d = (a1*dot(cross(p2_Hat,R1),p1_Hat) - dot(cross(p2_Hat,R2),p1_Hat) + a3*dot(cross(p2_Hat,R3),p1_Hat))/(a3*dot(cross(p2_Hat,p3_Hat),p1_Hat))    
    print "pd =", p1d,p2d,p3d

    p1 = Vector(p1d,RA1,Dec1)
    p2 = Vector(p2d,RA2,Dec2)
    p3 = Vector(p3d,RA3,Dec3)

    r1 = [p1[0]-R1[0],p1[1]-R1[1],p1[2]-R1[2]]
    r2 = [p2[0]-R2[0],p2[1]-R2[1],p2[2]-R2[2]]
    r3 = [p3[0]-R3[0],p3[1]-R3[1],p3[2]-R3[2]]
    rd = mag(r2)
    
    print "r =", r2

     # Correct for light speed change
    c = 173.1446 # speed of light in AU per day(JD)
    t1 = t1 - p1d/c
    t2 = t2 - p2d/c
    t3 = t3 - p3d/c

    # Recalculate tau after light speed change
    tau1 = k*(t1-t2)
    tau2 = k*(t3-t1)
    tau3 = k*(t3-t2)
    print "tau =", tau1,tau2,tau3

    # Recalculate R1, R2, and R3
    R1 = [R1[0]-(R1dot[0]*p1d)/c,R1[1]-(R1dot[1]*p1d)/c,R1[2]-(R1dot[2]*p1d)/c]
    R2 = [R2[0]-(R2dot[0]*p2d)/c,R2[1]-(R2dot[1]*p2d)/c,R2[2]-(R2dot[2]*p2d)/c]
    R3 = [R3[0]-(R3dot[0]*p3d)/c,R3[1]-(R3dot[1]*p3d)/c,R3[2]-(R3dot[2]*p3d)/c]
    print "R1 =", R1
    print "R2 =", R2
    print "R3 =", R3

    rdotB = VectorRDot(r1,r2,r3,rd,rdot,tau1,tau3)
    print "rdot =", rdotB

   


    f1 = f(rd,r2,rdot,tau1)
    g1 = g(rd,rdot,tau1)
    f3 = f(rd,r2,rdot,tau3)
    g3 = g(rd,rdot,tau3)
    
    print "f1 =",f1
    print "g1 =",g1
    print "f3 =",f3
    print "g3 =",g3
    
    
    # Initialize rdotA
    rdotA = [0,0,0]
    n = 1


    
    while mag([rdotB[0]-rdotA[0],rdotB[1]-rdotA[1],rdotB[2]-rdotA[2]]) > 10**-11: 

        # One iteration
        a1,a3 = a(rd,r2,rdotB,tau1,tau3)

        print " "
        print " "
        print "a1,a3 =", a1,a3

        p1d = (a1*dot(cross(R1,p2_Hat),p3_Hat) - dot(cross(R2,p2_Hat),p3_Hat) + a3*dot(cross(R3,p2_Hat),p3_Hat))/(a1*dot(cross(p1_Hat,p2_Hat),p3_Hat))
        p2d = (a1*dot(cross(p1_Hat,R1),p3_Hat) - dot(cross(p1_Hat,R2),p3_Hat) + a3*dot(cross(p1_Hat,R3),p3_Hat))/(-1*dot(cross(p1_Hat,p2_Hat),p3_Hat))
        p3d = (a1*dot(cross(p2_Hat,R1),p1_Hat) - dot(cross(p2_Hat,R2),p1_Hat) + a3*dot(cross(p2_Hat,R3),p1_Hat))/(a3*dot(cross(p2_Hat,p3_Hat),p1_Hat))    
        print "pd =", p1d,p2d,p3d

        p1 = Vector(p1d,RA1,Dec1)
        p2 = Vector(p2d,RA2,Dec2)
        p3 = Vector(p3d,RA3,Dec3)

        r1 = [p1[0]-R1[0],p1[1]-R1[1],p1[2]-R1[2]]
        r2 = [p2[0]-R2[0],p2[1]-R2[1],p2[2]-R2[2]]
        r3 = [p3[0]-R3[0],p3[1]-R3[1],p3[2]-R3[2]]
        rd = mag(r2)
        
        

##        c = 173.1446
##        t1 = t1 - p1d/c
##        t2 = t2 - p2d/c
##        t3 = t3 - p3d/c
##
##        tau1 = k*(t1-t2)
##        tau2 = k*(t3-t1)
##        tau3 = k*(t3-t2)
##        print "tau =", tau1,tau2,tau3

        print "R1 =", R1
        print "R2 =", R2
        print "R3 =", R3


        rdotA = VectorRDot(r1,r2,r3,rd,rdotB,tau1,tau3)
        print "r =", r2
        print "rdot =", rdotA

        f1 = f(rd,r2,rdot,tau1)
        g1 = g(rd,rdot,tau1)
        f3 = f(rd,r2,rdot,tau3)
        g3 = g(rd,rdot,tau3)
        
        print "f1 =",f1
        print "g1 =",g1
        print "f3 =",f3
        print "g3 =",g3
        print " "
        print " "



        # Next iteration
        a1,a3 = a(rd,r2,rdotB,tau1,tau3)
        
        print "a1,a3 =", a1,a3

        p1d = (a1*dot(cross(R1,p2_Hat),p3_Hat) - dot(cross(R2,p2_Hat),p3_Hat) + a3*dot(cross(R3,p2_Hat),p3_Hat))/(a1*dot(cross(p1_Hat,p2_Hat),p3_Hat))
        p2d = (a1*dot(cross(p1_Hat,R1),p3_Hat) - dot(cross(p1_Hat,R2),p3_Hat) + a3*dot(cross(p1_Hat,R3),p3_Hat))/(-1*dot(cross(p1_Hat,p2_Hat),p3_Hat))
        p3d = (a1*dot(cross(p2_Hat,R1),p1_Hat) - dot(cross(p2_Hat,R2),p1_Hat) + a3*dot(cross(p2_Hat,R3),p1_Hat))/(a3*dot(cross(p2_Hat,p3_Hat),p1_Hat))    
        print "pd =", p1d,p2d,p3d

        p1 = Vector(p1d,RA1,Dec1)
        p2 = Vector(p2d,RA2,Dec2)
        p3 = Vector(p3d,RA3,Dec3)

        r1 = [p1[0]-R1[0],p1[1]-R1[1],p1[2]-R1[2]]
        r2 = [p2[0]-R2[0],p2[1]-R2[1],p2[2]-R2[2]]
        r3 = [p3[0]-R3[0],p3[1]-R3[1],p3[2]-R3[2]]
        rd = mag(r2)

##        c = 173.1446
##        t1 = t1 - p1d/c
##        t2 = t2 - p2d/c
##        t3 = t3 - p3d/c
##
##        tau1 = k*(t1-t2)
##        tau2 = k*(t3-t1)
##        tau3 = k*(t3-t2)
##        print "tau =", tau1,tau2,tau3


        rdotB = VectorRDot(r1,r2,r3,rd,rdotA,tau1,tau3)
        print "r =", r2
        print "rdot =", rdotB


        f1 = f(rd,r2,rdotB,tau1)
        g1 = g(rd,rdotB,tau1)
        f3 = f(rd,r2,rdotB,tau3)
        g3 = g(rd,rdotB,tau3)  
        print "f1 =",f1
        print "g1 =",g1
        print "f3 =",f3
        print "g3 =",g3


    # Transform r and rdot into ecliptic coordinates
    eps = (23.5*pi)/180
    rmatrix = matrix([[r2[0]],[r2[1]],[r2[2]]])
    print rmatrix
    rdotmatrix = matrix([[rdotB[0]],[rdotB[1]],[rdotB[2]]])
    print rdotmatrix
    transformmatrix = matrix([[1,0,0],[0,cos(eps),sin(eps)],[0,-sin(eps),cos(eps)]])
    print transformmatrix

    rmatrix = transformmatrix*rmatrix
    print rmatrix
    rdotmatrix = transformmatrix*rdotmatrix
    print rdotmatrix

    r = [rmatrix[0][0],rmatrix[1][0],rmatrix[2][0]]
    rdot = [rdotmatrix[0][0],rdotmatrix[1][0],rdotmatrix[2][0]]
    print "r now equals", r
    print "rdot now equals", rdot

    return r, rdot


####################################################################################################

def Orbit_Determination(r,rdot):
    print " "
    print " "
    print "ORBITAL ELEMENTS"
    print " "

    # Semi-Major Axis, a
    k = 0.01720209895
    mu = 1
    a = 1/(2/mag(r) - (mag(rdot)**2)/mu)
    print "The semi-major axis equals", a

    # Eccentricity of Orbit, e
    e = sqrt(1 - (mag(cross(r,rdot))**2)/(mu*a))
    print "The eccentricity of the asteroid orbit equals", e

    # Inclination of Asteroid Orbit (i)
    l = cross(r,rdot)
    irads = atan(sqrt(l[0]**2 + l[1]**2)/l[2])
    i = irads*180/pi
    print "The inclination of the asteroid orbit equals", i

    # Longitude of Ascending Node (o)
    o1 = asin(l[0]/(mag(l)*sin(irads)))
    o2 = acos(-l[1]/(mag(l)*sin(irads)))
    orads = AngleAmb(o1,o2)
    o = orads*180/pi
    print "The longitude of the ascending node of the asteroid orbit equals", o

    # Argument of the Perihilion (w)
    U1 = asin(r[2]/(mag(r)*sin(irads)))
    U2 = acos((r[0]*cos(orads)+r[1]*sin(orads))/mag(r))
    U = AngleAmb(U1,U2)
    
    v1 = asin((((a*(1-e**2))/mag(l))*(dot(r,rdot)/mag(r)))/e)
    v2 = acos(((a*(1-e**2))/mag(r) - 1)/e)
    v = AngleAmb(v1,v2)
    
    wrads = U - v
    w = wrads*180/pi % 360
    print "The argument of the perihilion equals", w
    
    # Mean Anomaly (M)
    #n = sqrt(mu/(a**3))
    E = acos((1-mag(r)/a)/e)
    Mrads = E - e*sin(E)
    M = Mrads*180/pi
    print "The mean anomaly equals", M

    
    return a,e,i,o,w,M

###################################################################################################

# Input file as 7x9 array
fileString = raw_input("Please enter the name of the file: ")
pp = csv.reader(open(fileString,'rb'), delimiter = ',', quotechar = '|')
data = [[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]] # Create a 7x9 array
rowNum = 0 
for row in pp:
    colNum = 0
    for cell in row: 
        data[rowNum][colNum] = float(cell) 
        colNum += 1
    rowNum += 1
# Read in information from file
RA1 = data[1][0] + data[1][1]/60 + data[1][2]/3600
RA2 = data[4][0] + data[4][1]/60 + data[4][2]/3600
RA3 = data[7][0] + data[7][1]/60 + data[7][2]/3600
if data[1][3] > 0 or data[1][3] == 0:
    Dec1 = data[1][3] + data[1][4]/60 + data[1][5]/3600   
else:
    Dec1 = data[1][3] - data[1][4]/60 - data[1][5]/3600
       
if data[4][3] > 0 or data[1][3] == 0:
    Dec2 = data[4][3] + data[4][4]/60 + data[4][5]/3600
else:
    Dec2 = data[4][3] - data[4][4]/60 - data[4][5]/3600

if data[4][3] > 0 or data[4][3] == 0:
    Dec3 = data[7][3] + data[7][4]/60 + data[7][5]/3600
else:
    Dec3 = data[7][3] - data[7][4]/60 - data[7][5]/3600
    
R1 = [data[2][0],data[2][1],data[2][2]]
R2 = [data[5][0],data[5][1],data[5][2]]
R3 = [data[8][0],data[8][1],data[8][2]]
R1dot = [data[2][3],data[2][4],data[2][5]]
R2dot = [data[5][3],data[5][4],data[5][5]]
R3dot = [data[8][3],data[8][4],data[8][5]]
print RA1,RA2,RA3,Dec1,Dec2,Dec3

# Calculate time
year1 = data[0][2] 
year2 = data[3][2]
year3 = data[6][2]
month1 = data[0][1]
month2 = data[3][1]
month3 = data[6][1]
day1 = data[0][0]
day2 = data[3][0]
day3 = data[6][0]
ut1 = data[0][3] + data[0][4]/60 + data[0][5]/3600
ut2 = data[3][3] + data[3][4]/60 + data[3][5]/3600
ut3 = data[6][3] + data[6][4]/60 + data[6][5]/3600

# Convert time to Julian date
t1 = JulianDate(month1,day1,year1,ut1)
t2 = JulianDate(month2,day2,year2,ut2)
t3 = JulianDate(month3,day3,year3,ut3)
print "t1,t2,t3 =", t1,t2,t3
# Find r and rdot
r,rdot = Recursion(RA1,Dec1,t1,RA2,Dec2,t2,RA3,Dec3,t3,R1,R2,R3,R1dot,R2dot,R3dot)
print "r =", r
print "rdot =", rdot

# Find orbital elements
Orbit_Determination(r,rdot)
    
    
    
    


















