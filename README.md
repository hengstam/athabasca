# athabasca
 
/scripts:

	processRadarData.m:

		I haven't used this one.

	digitizeRadargram.m: 

		I didn't write this one. Run it, use it, and it'll output .mat files in /digitizedRadar that can be fed into the next script, cleanDigitizedRadar.m.

	cleanDigitizedRadar.m:

		You'll have to change the fileName variable on line 5 to process different files. Otherwise runs happily and will remove garbage values from a raw set of fresh picks. Outputs .mat files in the /scripts directory that the next script uses. Bonus: it displays what the new data looks like.

	generateVisdata.m:

		Produces the visdata.mat which the following script visualizes.

	visualizeAthabasca.m: 

		Visualizes the whole glacier based on the visdata.mat file produced by the previous script.

visdata.mat:
    
    longs (degrees), lats (degrees), elevs (m), times (matlab time), sources (see below):
        These fields are the cleaned and linearly interpolated data obtained by the handheld GPS. The interpolation was performed to eliminate the "bursty" sort of data collection the raw data displayed. Interpolation is performed by generating an artifical and evenly-spaced time series, then interpolating existing position data to obtain position data that would correspond to the artificial times. 
        All further position data have the same units. These data and all further position data are sorted such that time is monotonically ordered.

    longPicks, latPicks, elevPicks, timePicks, sourcePicks: 
        These fields are the data of the surface position of each bedrock pick, based on the interpolated handheld GPS data. Times are found by assuming that each shot location corresponds to one of the artifical times developed in interpolating the above data. Where a bedrock pick falls between two shots, the time of that bedrock pick is adjusted to reflect this intermediate value.

    longsCSRS, latsCSRS, elevsCSRS, timesCSRS:
        These are the the raw data obtained by the CSRS GPS.

    longCSRSPicks, latCSRSPicks, elevCSRSPicks, timeCSRSPicks, sourceCSRSPicks:
        The times contained in timePicks are used to generate these surface positions of bedrock picks based on the CSRS data. Where the bedrock time picks fall between CSRS data collection times, linear interpolation is used to generate intermediate position data.

    tri, triPick, triCSRS, triCSRSPick:
        These are Delauny triangulation meshes (https://www.mathworks.com/help/matlab/ref/delaunay.html), representing the set of triangles that make up the triangulation. Required to visualize the Delauny meshes.
        Each matrix is a triangulation mesh for the corresponding datasets of lats and longs.

    zPicks: 
        These are the raw depths (in meters) of the bedrock picks.

    sources, sourcePicks:
        Indicate the file origin of the datapoint. Integer values correspond to the character names stored in the cell array `sourceNames`.

    sourceNames:
        Serves as a legend for `sources` and `sourcePicks`.