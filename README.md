# athabasca
 
/scripts:

	processRadarData.m:

		I haven't used this one.

	digitizeRadargram.m: 

		I didn't write this one. Run it, use it, and it'll output .mat files in /digitizedRadar that can be fed into the next script, cleanDigitizedRadar.m.

	cleanDigitizedRadar.m:

		You'll have to change the fileName variable on line 5 to process different files. Otherwise runs happily and will remove garbage values from a raw set of fresh picks. Outputs .mat files in the /scripts directory that the next script uses. Bonus: it displays what the new data looks like.

	generateVisdata.m:

		Produces the datafile which the following script visualizes.

	visualizeAthabasca.m: 

		Visualizes the whole glacier. Relies on the previous script to run.