This set of programs bundled as uwsDropletEvaporationCode can be used to 
predict the diameter, temperature and mass fraction histories of an 
evaporating urea-water-solution droplet. 

Copyright (C) 2021 Mikhil Surendran
     
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

The method used in uwsDropletEvaporationCode is mostly based on the 
single component droplet evaporation model proposed by Abramzon and 
Sirignano (1989), and the details of this method have been published and 
are available online.
(S. Mikhil, S. Bakshi, T. N. C. Anand, A computational model for the 
evaporation of urea-water-solution droplets exposed to a hot air stream, 
International Journal of Heat and Mass Transfer 168 (2021) 120878.  
10.1016/j.ijheatmasstransfer.2020.120878).

The code consists of a main function, an ode45 function, an input file 
and several other functions to evaluate the required thermophysical 
properties, non-dimensional numbers etc. In order to use the program for 
predicting the evaporation behavior of UWS droplets, or other 
single-component droplets, the values of the variables given in the 
input file should be selected appropriately and the main function should 
be compiled. The program can be readily used to predict the evaporation 
behavior of any single-component liquid droplet by incorporating the 
required property values of the desired liquids, in the respective 
functions, and by setting the mass fraction of component two to zero. 
The program may also be used for other bi-component liquids by making 
minor modifications to the code made available here. 
