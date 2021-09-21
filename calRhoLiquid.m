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

%% This function calculates the density of individual components in the liquid.
%  Film temperature is evaluated using the weighting parameter defined as the "FilmRulePara"
%  FilmRulePara is 0 for liquids. 
%  Inputs and outputs are in SI units

function [density, dRhoBydTemperature] = calRhoLiquid(ambientPressure, ambientTemperature, initialDropletTemperature, boilingPoint, surfaceTemperature, filmRulePara, liquid, salt, saltMassFraction, calMethod, propCharacter, expOrCorrln)

%% select the liquid
water = strcmp(liquid,'water');
moltenUrea = strcmp(liquid,'urea');
nHeptane = strcmp(liquid,'nHeptane');
nDecane = strcmp(liquid,'nDecane');

%% check if it is urea water solution 
urea = strcmp(salt,'urea') ;

%% check if the density of UWS is to be calculated based on the density of solid urea 
moltenUreaDensity = strcmp(expOrCorrln,'moltenUrea');
constant = strcmp(expOrCorrln,'constant');  % density of solid urea is taken through out and is independent of temperature. 

%% select the correlation
mikhil = strcmp(calMethod, 'mikhil');

%% dependence on temperature 
temperatureDependent = strcmp (propCharacter, 'temperatureDependent') ; 
temperatureIndependent = strcmp (propCharacter, 'temperatureIndependent') ;

%% set the temperature at which properties are evaluated 
if temperatureDependent == 1
    
    if le(surfaceTemperature, boilingPoint)
        
        temperature = surfaceTemperature ; 
        
    else
        
        temperature = boilingPoint ;    % if the evaluated droplet temperature is higher than the boiling point of a component, then its properties are evaluated at its boiling point. 
        
    end
    
    if moltenUrea == 1                        % density of urea is evaluated at the actual temperature
        
        temperature = surfaceTemperature ; 
        
    end
    
end

if temperatureIndependent == 1 
    
    temperature = 0.5 * (initialDropletTemperature + boilingPoint) ;   %  Abramzon evaluated liquid properties at an average temperature (To + Tboil)/2  
    
end

%% calculate the liquid density 

if mikhil == 1  
    
    if moltenUrea == 1 && constant == 1
        
        density = 1335 ;
        
        dRhoBydTemperature = 0 ; 

% -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------                
    elseif moltenUrea == 1 && moltenUreaDensity == 1 && temperature < 406.15  % this subroutine uses the density of solid urea to evaluate the density of UWS below the melting point of urea. Above the melting point it uses a correlation for molten urea. 
        
        density = 1335 ;
        
        dRhoBydTemperature = 0 ;   
        
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------                
    else   % these subroutines evaluate the liquid phase densities of the components
    
        if water == 1
   
            A = 0.34710;
            B = 0.274;
            n = 0.28571;
            criticalTemperature = getCriticalTemperature('water'); 
    
        elseif moltenUrea == 1
        
            A = 0.56982;
            B = 0.337;
            n = 0.28571;
            criticalTemperature = getCriticalTemperature('urea');
        
        elseif nHeptane == 1
        
            A = 0.23237;
            B = 0.2602;
            n = 0.2791;
            criticalTemperature = getCriticalTemperature('nHeptane');
    
        elseif nDecane == 1

            A = 0.23276;
            B = 0.2524;
            n = 0.2857;
            criticalTemperature = getCriticalTemperature('nDecane');
    
        end   
    
        density = A * B ^(-(1-(temperature/criticalTemperature))^n) * 1000 ; % kg/m^3   % correlations used by abramzon wasn't available so used the correlation given in Chemical Properties Handbook by Carl Yaws
    
        if temperatureIndependent == 1 

            dRhoBydTemperature = 0 ; 
        
        else

            dRhoBydTemperature = (n * log(B) / criticalTemperature) * (1 - (temperature/criticalTemperature))^(n-1) * A * B ^(-(1-(temperature/criticalTemperature))^n) * 1000 ; % derivative of density with respect to temperature
    
        end
        
    end
    
end

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
