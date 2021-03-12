# astro-tools

A collection of functions as I collect useful equations and algorithms for astrodynamics and space research! This will probably eventually become several repositories eventually, but for the moment, having a single file with functions in it for different applications seems sufficient.

Current tools are for 
 - The two-body problem (twobodytools.py)
 - The CR3BP (threebodytools.py)
 - Spacecraft attitude dynamics (attitudedynamics.py)
 - Planetary constants (bodies.py)
 
To use these functions, install them somewhere on your PYTHONPATH (or add the location of the files to the PYTHONPATH) (or put them in the same directory where the script from which you're calling the functions is in) and use `import functionfile` or variants thereof within your Python script
