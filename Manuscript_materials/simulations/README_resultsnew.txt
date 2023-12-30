Results from long run; ./results_new/

We get 5 35-years long slices (0-35y, 50-85y, 250-285y, 1000-1035y, and 2465-2500y)

binaout_p_S_Y_id: kept the occupancy pattern of year Y, for period p, species S, in simulation with id = id; binary matrix (nrep x nsites)
bysite_p_S_Y_id: give the average occupancy by site , for period p, species S, slice Y and id = id; real matrix (nrep x nsites) 


Init = Spring 83										Full Matrix		Interactions 			172507
Init = Spring 83								 		Full Matrix	 	Interactions 			316018	+ offset (0.0005)
Init = Spring 83										Full Matrix		No interactions 		446512
Init = Spring 83								 		Full Matrix	 	No interactions 		337776	+ offset (0.0005)

Init = (Summ 82 + Spring 83 - Summ 82 * Spring 83) 		Full Matrix 	Interactions 			656864
Init = (Summ 82 + Spring 83 - Summ 82 * Spring 83)		Full Matrix 	No interactions 		635880

Init = Random 0.2								 		Full Matrix 	Interactions 			660193
Init = Random 0.2								 		Full Matrix	 	Interactions 			150370	+ offset (0.0005)
Init = Random 0.2										Full Matrix 	No interactions 		405678
Init = Random 0.2								 		Full Matrix	 	No interactions 		278148	+ offset (0.0005)


Init = Spring 83										No space		Interactions 			205143
Init = Spring 83								 		No space	 	Interactions 			844859 + offset (0.0005)
Init = Spring 83										No space		No interactions 		819821
Init = Spring 83										No space		No interactions 		850617 + offset (0.0005)

Init = (Summ 82 + Spring 83 - Summ 82 * Spring 83) 		No space	 	Interactions 			13644
Init = (Summ 82 + Spring 83 - Summ 82 * Spring 83)		No space		No interactions 		74757

Init = Random 0.2								 		No space		Interactions 			837704
Init = Random 0.2								 		No space	 	Interactions 			410908 + offset (0.0005)
Init = Random 0.2										No space	 	No interactions 		703710
Init = Random 0.2								 		No space	 	No interactions 		526999 + offset (0.0005)

########################################################################################################################

Files : 

	- 	allSummResults.RDS : 		All results are concatenate in nested list ; contain both all replication for all situation (dt list), and averag long term prob. (dt_longterm)
									It does not contain the bianout_x_x_x_x entries however. 

	- 	allresults.zip : 			Compressed results files. 


########################################################################################################################

Plots: 

i) Matrice complète de distance + initialisation comme en 83 + interactions + summer period + longterm (2500 ans)
ii) Matrice complète de distance + initialisation random (=0.2) + interactions + summer period + longterm (2500 ans)

iii) Matrice complète de distance + initialisation comme en 83 + interactions + summer period + "longterm" (35 ans)
iv) Matrice complète de distance + initialisation random (=0.2) + interactions + summer period + "longterm" (35 ans)

v) Pas d'espace explicit + initialisation comme en 83 + interactions + summer period + "longterm" (2500 ans)
vi) Pas d'espace explicit + initialisation random (=0.2) + interactions + summer period + "longterm" (2500 ans) 

Versions (XX-b): Without interactions.
Version (xx)-NoE: Removing extinct MP