% This set of programs bundled as uwsDropletEvaporationCode can be used to 
% predict the diameter, temperature and mass fraction histories of an 
% evaporating urea-water-solution droplet. 

% Copyright (C) 2021  Mikhil Surendran
     
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% email: mikhilsuren@gmail.com

%% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% The method used in uwsDropletEvaporationCode is mostly based on the single 
% component droplet evaporation model proposed by Abramzon and Sirignano (1989), 
% and the details of this method have been published and are available 
% online. (Mikhil et al. 10.1016/j.ijheatmasstransfer.2020.120878).

% The code consists of a main function, an ode45 function, an input file 
% and several other functions to evaluate the required thermophysical 
% properties, non-dimensional numbers etc. The program can also be readily 
% used to predict the evaporation behavior of any single component liquid 
% droplet by incorporating the required property values of the desired 
% liquids, and may be used for other bi-component liquids by making 
% minor modifications to the code. 

%% This function defines the critical pressures. 
%  Critical temperature of the liquid is used in the evalution of certain 
%  properties like density of the liquid phase.

function [criticalPressure] = getCriticalPressure(liquid)

nitrogen = strcmp(liquid, 'N2');
air = strcmp(liquid, 'air');
water = strcmp(liquid, 'water');
urea = strcmp(liquid, 'urea');
nHeptane = strcmp(liquid, 'nHeptane');
nDecane = strcmp(liquid, 'nDecane');
waterVapor = strcmp(liquid, 'waterVapor');
ureaVapor = strcmp(liquid, 'ureaVapor');
nHeptaneVapor = strcmp(liquid, 'nHeptaneVapor');
nDecaneVapor = strcmp(liquid, 'nDecaneVapor');

if nitrogen == 1
    
    criticalPressure = 34 * 10^5 ; 
    
elseif air == 1
    
    criticalPressure = 37.9 * 10^5 ; 

elseif water == 1 || waterVapor == 1
    
    criticalPressure = 220.6 * 10^5 ;
    
elseif nHeptane == 1 || nHeptaneVapor == 1
    
    criticalPressure = 27.35 * 10^5 ;
   
elseif nDecane == 1 || nDecaneVapor == 1
    
    criticalPressure = 21.05 * 10^5 ;  
    
end
 