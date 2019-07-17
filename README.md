# OrbitDetermination
This software uses Gauss’s Method of Orbital Determination to determine the six orbital elements of an asteroid, given the data from three separate observations made by the same telescope at different dates. Originally made for SSP 2012, Westmont site.

INPUT: a csv file containing a table that holds, for each Observation, the Right Ascension (hours, minutes, seconds / degrees) and the Declination (degrees, minutes, seconds / degrees) of the asteroid observed, as well as the time (UT / Julian Date) of the observation. See data.csv for example formatting.

OUTPUT: 
	1. The Semi-Major Axis (a), in AU
	2. The Eccentricity (e)
	3. The Inclination (i), in degrees
	4. The Longitude of the Ascending Node (capital omega), in degrees
	5. The Argument of the Periapsis, (omega), in degrees
	6. The Mean Anomaly (M_0), in degrees

To run the Orbit Determination program with the example data (which comes from Observations of asteroid 1998 QS52 made 6/18/12, 7/15/12, and 7/19/12 at Westmont College using the Mead 14" and Keck 12" telescopes), run OrbitDetermination.py using VPython in IDLE 2.7.9 and input "data.csv" when prompted for the input file.

This program was written by E.J. Khouri in 2012.
