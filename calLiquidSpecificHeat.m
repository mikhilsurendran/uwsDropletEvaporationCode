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

%% This function calculates the specific heat of the liquid phase. It may be suitably modified to use other expressions. 
%  Film temperature is evaluated using the weighting parameter defined as the "FilmRulePara"
%  FilmRulePara = 0 for liquids

function [specificHeat] = calLiquidSpecificHeat(ambientPressure, ambientTemperature, initialDropletTemperature, boilingPoint, surfaceTemperature, filmRulePara, liquid, calMethod, propCharacter)

%% select the liquid
water = strcmp(liquid,'water');
urea = strcmp(liquid,'urea');
nHeptane = strcmp(liquid,'nHeptane');
nDecane = strcmp(liquid,'nDecane');

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
    
    if urea == 1                        % specific heat of urea is evaluated at the actual temperature
        
        temperature = surfaceTemperature ; 
        
    end
    
end

if temperatureIndependent == 1 
    
    temperature = 0.5 * (initialDropletTemperature + boilingPoint) ;  % as done in Abramzon's paper
    
end

%% Calculate specific heat of the liquid

if mikhil == 1
    
    if water == 1     % correlations given in Chemical Properties Handbook by Carl L. Yaws
    
        A = 92.053;
        B = -3.9953*10^(-2);
        C = -2.1103*10^(-4);
        D = 5.3469*10^(-7);
    
        molecularWeight = 18.02; % g/mol
    
    elseif urea == 1
    
        A = 965.507;
        B = -5.0993;
        C = 1.0028*10^(-2);
        D = -6.3799*10^(-6); 
    
        molecularWeight = 60.06; % g/mol
    
    elseif nHeptane == 1
    
        A = 101.121;
        B = 9.7739 * 10^(-1);
        C = -3.0712 * 10^(-3);
        D = 4.1844 * 10^(-6); 
    
        molecularWeight = 100.204; % g/mol
    
    elseif nDecane == 1
    
        A = 79.741;
        B = 1.6926;
        C = -4.5287 * 10^(-3);
        D = 4.9769 * 10^(-6); 
    
        molecularWeight = 142.285; % g/mol
    
    end
  
        specificHeat = A + B * temperature + C * temperature^2 + D * temperature^3;   %J/mol K 
    
        specificHeat = specificHeat / (molecularWeight / 1000); % J/kgK   
    
end