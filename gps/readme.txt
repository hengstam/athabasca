readme for GPS folder
03 june 2019

radarCoordinatesAndTime_handheldGPS.csv - coordinates for the heldheld GPS moving with the radar sleds.
First column = shot number, a unique number for each depth sounding.
Second col = transect number, a unique number for each cross-glacier or down-glacier profile
Third col = latitude (degree)
Fourth col = longitude (degree)
Fifth col = height (m above sea level)

raw_gps_data - you will likely not have to use this folder that contains raw gps data (not in easily read formats)

csrs - processed GPS data output from the National Resources Canada's Canadian Spatial Reference System (CSRS).
This contains the files:
ath<xxxx>.pdf – summary file showing plots of satellite positions, antenna coordinates, and data quality during the transect.
ath<xxxx>.csv – a text file containing antenna coordinates. This will be your main file to work with. There is a header that explains what each column represents.
ath<xxxx>.sum – summary text file describing the survey. Probably not very useful.

File naming convention: ath<DOY><SURVEY_NUM> where DOY is the day of year (where DOY 1 = Jan 1, DOY 365 = Dec 31); SURVEY_NUM is the sequential survey number that was run on that DOY.

More on CSRS at:
https://www.nrcan.gc.ca/earth-sciences/geomatics/geodetic-reference-systems/tools-applications/10925

 