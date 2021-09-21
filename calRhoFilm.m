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

%% This function calculates the density of the mixture in the film near the droplet surface using ideal gas law. It may be suitably modified to use other laws. 

function [density] = calRhoFilm(ambientPressure, ambientTemperature, surfaceTemperature, molecularWeightOne, molecularWeightTwo, vaporOneMoleFraction, vaporTwoMoleFraction, ambGasMoleFraction, vapOneDensity, vapTwoDensity, ambGasDensity, vapOneMassFraction_film, vapTwoMassFraction_film, gasMassFraction_film, filmRulePara, gas, rule, calMethod)

% select the correlation
mikhil = strcmp(calMethod, 'mikhil');

%% Molar mass of the ambient gas (nitrogen or air)
nitrogen = strcmp(gas,'N2');
air = strcmp(gas,'air');

if nitrogen == 1
    
    gasMolecularWeight = 28.01 ; % g/mol
    
elseif air == 1
    
    gasMolecularWeight = 28.97 ; % g/mol
    
end

%% Calculate the mixture density

if mikhil == 1 

    % Molar mass of the mixture
    molarMassMixture = vaporOneMoleFraction*molecularWeightOne + vaporTwoMoleFraction*molecularWeightTwo + ambGasMoleFraction*gasMolecularWeight ; 
    

    %  density is calculated based on the selected gas law
    ideal = strcmp(rule,'ideal');
    real = strcmp(rule,'real');

    if ideal == 1
    
        filmTemperature = surfaceTemperature*(1-filmRulePara) + ambientTemperature*filmRulePara;
    
        gasConstant = 8.314 ; % J/K mol
    
        density = ambientPressure*molarMassMixture/(1000*gasConstant*filmTemperature); % kg/m^3
    
    end

    if real == 1
        
        error('Density evaluation for vapor film not defined when using real gas laws');
        
    end

end



