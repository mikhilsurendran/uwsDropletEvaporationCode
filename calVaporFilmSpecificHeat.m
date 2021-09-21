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

%% This function calculates the specific heat of the vapour film surrounding the droplet.  
%  Another function calVaporSpecificHeat(...) is used to calculate
%  individual specific heats of the vapors/gasses constituting the film.

function [specificHeat] = calVaporFilmSpecificHeat(specificHeatOne, specificHeatTwo, gasSpecificHeat, compOneMassFraction_surf, compOneMassFraction_inf, compTwoMassFraction_surf, compTwoMassFraction_inf, filmRulePara, calMethod)

% select the correlation
mikhil = strcmp(calMethod, 'mikhil');

%% calculate the film specific heats

if mikhil ==1
    
    compOneMassFraction_film = compOneMassFraction_surf*(1-filmRulePara) + compOneMassFraction_inf*filmRulePara ;

    compTwoMassFraction_film = compTwoMassFraction_surf*(1-filmRulePara) + compTwoMassFraction_inf*filmRulePara ;

    gasMassFraction_film = 1 - (compOneMassFraction_film + compTwoMassFraction_film) ; 
    
    specificHeat = compOneMassFraction_film*specificHeatOne + compTwoMassFraction_film*specificHeatTwo + gasMassFraction_film*gasSpecificHeat ;  % J/kg K
    
end

