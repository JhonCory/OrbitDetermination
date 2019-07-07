from __future__ import division
from math import *
from visual import *
import numpy
import csv

###############################################################################

"""
Orbit Determination
Written by John Khouri

Initial revision: 20 July 2012 (John Khouri)
Last revision: 6 July 2019 (John Khouri)

The purpose of this program is to use the Method of Gauss to generate the six
orbital elements of an asteroid using the data from three separate
observations of this asteroid.

The inputs include (for each observation) the time, the right ascension and
declination of the asteroid, and the positional and velocity vectors for the
sun relative to the earth at the time of the observation.
"""

###############################################################################
###############################################################################

# Convert degrees to radians
def degs2rads(degs):
    return (degs/360)*2*pi
    
# Convert hours to radians
def hrs2rads(hrs):
    return (hrs*pi)/12
    
# Take the dot product of two vectors
def dot(a,b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

# Take the cross product of two vectors
def cross(a,b):
    return [a[1]*b[2]-a[2]*b[1],-a[0]*b[2]+a[2]*b[0],a[0]*b[1]-a[1]*b[0]]

# Return the sum of two vectors
def add(a,b):
    return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]

# Return the difference between two vectors
def subtract(a,b):
    return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]

# Apply a scalar to a vector
def apply_scalar(x, v):
    return [x*v[0], x*v[1], x*v[2]]

# Given the RA and Dec, return the corresponding unit vector
def get_unit_vector(RA,Dec):
    return [cos(RA)*cos(Dec), sin(RA)*cos(Dec), sin(Dec)]

# Calculate a vector from its magnitude and unit vector
def get_vector(mag, unitVector):
    return [mag*unitVector[0], mag*unitVector[1], mag*unitVector[2]]

# Calculate the distances from the asteroid for each observation
def get_distances(a1, a3, p1_hat, p2_hat, p3_hat, R1, R2, R3):
    p1d = (a1*dot(cross(R1,p2_hat),p3_hat) - dot(cross(R2,p2_hat),p3_hat) + a3*dot(cross(R3,p2_hat),p3_hat))/(a1*dot(cross(p1_hat,p2_hat),p3_hat))
    p2d = (a1*dot(cross(p1_hat,R1),p3_hat) - dot(cross(p1_hat,R2),p3_hat) + a3*dot(cross(p1_hat,R3),p3_hat))/(-1*dot(cross(p1_hat,p2_hat),p3_hat))
    p3d = (a1*dot(cross(p2_hat,R1),p1_hat) - dot(cross(p2_hat,R2),p1_hat) + a3*dot(cross(p2_hat,R3),p1_hat))/(a3*dot(cross(p2_hat,p3_hat),p1_hat))    
    return p1d, p2d, p3d

# Calculate the positional vectors of the asteroid relative to the Earth
# for each observation
def get_pos_vectors(p1d, p2d, p3d, p1_hat, p2_hat, p3_hat, R1, R2, R3):
    p1 = get_vector(p1d, p1_hat)
    p2 = get_vector(p2d, p2_hat)
    p3 = get_vector(p3d, p3_hat)

    r1 = subtract(p1, R1)
    r2 = subtract(p2, R2)
    r3 = subtract(p3, R3)
    return r1, r2, r3

def f(rd,r2,rdot,tau):
    return 1 - (tau**2)/(2*rd**3) + ((tau**3)*dot(r2,rdot))/(2*(rd**5)) + ((tau**4)/24)*((3/(rd**3))*((dot(rdot,rdot)/(rd**2))-1/(rd**3))-(15/(rd**2))*dot(rdot,rdot)**2+1/(rd**6))
    return f

def g(rd,rdot,tau):
    return tau - (tau**3)/(6*rd**3) + ((tau**4)*dot(rdot,rdot))/(4*(rd**5))

# Recalculate rdot, the velocity vector of the asteroid relative to the sun,
# using previous estimation
def calc_rdot(r1,r2,r3,rdot_old,tau1,tau3):
    rd = mag(r2)
    f3 = f(rd,r2,rdot_old,tau3)
    g3 = g(rd,rdot_old,tau3)
    r1dot = apply_scalar(1/g3, subtract(r3, apply_scalar(f3, r2)))
    f1 = f(rd,r2,rdot_old,tau1)
    g1 = g(rd,rdot_old,tau1)
    r2dot = apply_scalar(1/g1, subtract(r1, apply_scalar(f1, r2)))
    rdot_new = apply_scalar(1/2, add(r1dot, r2dot))
    return rdot_new, rdot_old

# Reapproximate a1 and a3 using the most recent estimations
def reapprox_a(r2,rdot,tau1,tau3):
    rd = mag(r2)
    common_divisor = \
	f(rd,r2,rdot,tau1) * g(rd,rdot,tau3) - \
	f(rd,r2,rdot,tau3) * g(rd,rdot,tau1)
    a1 = g(rd,rdot,tau3) / common_divisor
    a3 = -g(rd,rdot,tau1) / common_divisor
    return a1,a3

# Recalculate tau1, tau2, tau3 based on latest estimates for distance of
# the earth from the asteroid 
def correct_for_light(t1, t2, t3, p1d, p2d, p3d):
    # c is the speed of light, expressed in AU per day (JD)
    c = 173.1446

    # Factor in the speed of light and the current best estimates of the 
    # distances from the asteroid to recalculate observed times from observation
    # times
    #   
    # This is important because if the observation of an object took place at
    # time t, then the observed position would be from time t - d/c, where d is
    # the distance to the object and c is the speed of light
    obd_t1, obd_t2, obd_t3 = t1-p1d/c, t2-p2d/c, t3-p3d/c

    # Recalculate tau after light speed adjustment
    tau1, tau2, tau3 = get_tau(obd_t1, obd_t2, obd_t3)
    return tau1, tau2, tau3

def get_tau(t1, t2, t3):
    # k is the Gaussian gravitational constant, expressed in radians per day
    # NOTE - this program was written July 2012. At the time, the standard
    # value for the Gaussian gravitational constant was the one defined by
    # Carl Gauss in 1809. However, this may no longer be the acceptable value
    # for k. If not, the value below for k can be replaced.
    # TODO - make global variable?
    k = 0.01720209895

    tau1 = k*(t1-t2)
    tau2 = k*(t3-t1)
    tau3 = k*(t3-t2)
    return tau1, tau2, tau3

def julian_date(month,day,year,time):
    J0 = 367*year - int((7*(year + int((month+9)/12)))/4) + int((275*month)/9) + day + 1721013.5
    JD = J0 + time/24
    return JD

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
    
###############################################################################

"""
Recurse over three observational inputs using Gauss's Method and estimate the
position and velocity vectors of the asteroid relative to the earth at
the time associated with the second observation: r and rdot
"""
def calculate_vectors(RA1,Dec1,t1,RA2,Dec2,t2,RA3,Dec3,t3,R1,R2,R3,R1dot,R2dot,R3dot):
    # c is the speed of light, expressed in AU per day (JD)
    c = 173.1446

    # Convert RA and Dec into radians
    RA1, RA2, RA3 = hrs2rads(RA1), hrs2rads(RA2), hrs2rads(RA3)
    Dec1, Dec2, Dec3 = degs2rads(Dec1), degs2rads(Dec2), degs2rads(Dec3)

    # For each observation, calculate the unit vector of p, the position vector
    # of the asteroid relative to the earth
    p1_hat = get_unit_vector(RA1, Dec1)
    p2_hat = get_unit_vector(RA2, Dec2)
    p3_hat = get_unit_vector(RA3, Dec3)
    
    # Calculate initial values for scalars tau
    tau1, tau2, tau3 = get_tau(t1, t2, t3)

    # Make initial estimations of scalars a1 and a3
    a1 = tau3/tau2
    a3 = -tau1/tau2

    # For each vector p (position vector of the asteroid relative to the earth),
    # calculate its magnitude pd, i.e. the distance from the asteroid
    p1d, p2d, p3d = get_distances(a1, a3, p1_hat, p2_hat, p3_hat, R1, R2, R3)

    tau1, tau2, tau3 = correct_for_light(t1, t2, t3, p1d, p2d, p3d)

    # Now that decent estimations of p1d, p2d, p3d have been made, recalculate
    # R1, R2, and R3 based on light speed
    # This will only be done once so that the variance in values calculated
    # for rdot may converge to zero
    R1 = subtract(R1, apply_scalar(p1d/c, R1dot))
    R2 = subtract(R2, apply_scalar(p2d/c, R2dot))
    R3 = subtract(R3, apply_scalar(p3d/c, R3dot))

    r1, r2, r3 = get_pos_vectors(p1d, p2d, p3d, p1_hat, p2_hat, p3_hat, R1, R2, R3)
         
    # Use linear approximation to approximate vector rdot (the velocity vector
    # of the asteroid relative to the sun around the second observation)
    rdot_new = apply_scalar(1/(tau3-tau1), subtract(r3,r1))

    print "r (initial) = ", r2
    print "rdot (initial) = ", rdot_new
    print ""

    rdot_prev = [0,0,0]
    n = 1
    # Iterate until the values for rdot converge
    while mag([rdot_new[0]-rdot_prev[0],rdot_new[1]-rdot_prev[1],rdot_new[2]-rdot_prev[2]]) > 10**-11: 
        # Recalculate a1 and a3
        a1,a3 = reapprox_a(r2,rdot_new,tau1,tau3)

        p1d, p2d, p3d = get_distances(a1, a3, p1_hat, p2_hat, p3_hat, R1, R2, R3)
        r1, r2, r3 = get_pos_vectors(p1d, p2d, p3d, p1_hat, p2_hat, p3_hat, R1, R2, R3)

        tau1, tau2, tau3 = correct_for_light(t1, t2, t3, p1d, p2d, p3d)

        # Recalculate rdot
        rdot_new, rdot_prev = calc_rdot(r1,r2,r3,rdot_new,tau1,tau3)
        print "r = ", r2
        print "rdot = ", rdot_new
        print ""

    # Transform r and rdot into ecliptic coordinates
    eps = (23.5*pi)/180
    rmatrix = matrix([[r2[0]],[r2[1]],[r2[2]]])
    rdotmatrix = matrix([[rdot_new[0]],[rdot_new[1]],[rdot_new[2]]])
    transformmatrix = matrix([[1,0,0],[0,cos(eps),sin(eps)],[0,-sin(eps),cos(eps)]])

    rmatrix = transformmatrix*rmatrix
    rdotmatrix = transformmatrix*rdotmatrix

    r = [rmatrix[0][0],rmatrix[1][0],rmatrix[2][0]]
    rdot = [rdotmatrix[0][0],rdotmatrix[1][0],rdotmatrix[2][0]]
    return r, rdot


###############################################################################

"""
Using the position and velocity vectors of the asteroid relative
to the sun (r and rdot), calculate and print its six orbital
elements.
"""
def determine_orbital_elements(r,rdot):
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
    E = acos((1-mag(r)/a)/e)
    Mrads = E - e*sin(E)
    M = Mrads*180/pi
    print "The mean anomaly equals", M

###############################################################################

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

"""
Read in information from file

Consider rows 1-3, which contain the data for the first
of three observations. This is how the data is organized:

day      month     year      UT (hrs)  UT (mins)    UT (secs)
RA (hrs) RA (mins) RA (secs) Dec (deg) Dec (arcmin) Dec (arcsec)
R [i]    R [j]     R [k]     Rdot [i]  Rdot [j]     Rdot [k]

1. The first row contains the day/month/year/time (in UT) that
the observation was made.
2. The second row contains the right ascension and declination
of the asteroid at the time of the observation.
3. The third row contains the values (i, j, and k) for the vectors
"R" and "R dot," representing the sun's position and velocity
relative to the earth at the time of the observation.

"""


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
# Convert times to Julian dates
t1 = julian_date(month1,day1,year1,ut1)
t2 = julian_date(month2,day2,year2,ut2)
t3 = julian_date(month3,day3,year3,ut3)


# Find r and rdot
r,rdot = calculate_vectors(RA1,Dec1,t1,RA2,Dec2,t2,RA3,Dec3,t3,R1,R2,R3,R1dot,R2dot,R3dot)

# Find orbital elements
determine_orbital_elements(r,rdot)
