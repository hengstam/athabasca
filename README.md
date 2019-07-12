# athabasca
 
/scripts:

	processRadarData.m:

		I haven't used this one.

	digitizeRadargram.m: 

		I didn't write this one. Run it, use it, and it'll output .mat files in /digitizedRadar that can be fed into the next script, cleanDigitizedRadar.m.

	cleanDigitizedRadar.m:

		You'll have to change the fileName variable on line 5 to process different files. Otherwise runs happily and will remove garbage values from a raw set of fresh picks. Outputs .mat files in the /scripts directory that the next script uses. Bonus: it displays what the new data looks like.

	visualizeAthabasca.m: 

		No input require. Just run it and it'll display everything you've already processed with cleanDigitizedRadar, provided the output is in the same directory. Might take a bit to finish running.