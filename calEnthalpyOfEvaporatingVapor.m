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

%% This function calculates the enthalpy of vaporization of the evaporating vapor. 
%  Enthalpy of vaporization is evaluated as a weighted average of
%  individual enthalpies, using mole fractions as the weights. 
%  Another function calVaporizationEnthalpy(...) is used to calculate
%  individual vaporization enthalpies of the components in the liquid.

function [enthalpy] = calEnthalpyOfEvaporatingVapor(enthalpyOne, enthalpyTwo, compOneMoleFraction, compTwoMoleFraction, molecularWeightOne, molecularWeightTwo, calMethod)

%% select the correlation
mikhil = strcmp(calMethod, 'mikhil');

%% calculate the enthalpy of the evaporating mixture

if mikhil == 1 
    
    if (compOneMoleFraction > 0 || compTwoMoleFraction > 0 )
    
        weightOne = compOneMoleFraction / (compOneMoleFraction + compTwoMoleFraction);
        weightTwo = 1 - weightOne ;
    
    else 
    
        weightOne = 0;
        weightTwo = 0;
    
    end

    enthalpy = enthalpyOne * weightOne + enthalpyTwo * weightTwo; % J/mol

    vaporMolecularWeight = (molecularWeightOne * weightOne + molecularWeightTwo * weightTwo) / 1000 ;% kg/mol 

    if vaporMolecularWeight > 0
    
        enthalpy = enthalpy / vaporMolecularWeight; % J/kg
    
    else 
    
        enthalpy = 0; 
    
    end
    
end
