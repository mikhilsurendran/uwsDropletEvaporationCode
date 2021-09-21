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

%% This function calculates the boiling points of the components of the bi-component liquid / solution. 
%  This boiling point may then be used for evaluating other liquid phase properties (Abramzon and Sirignano evaluated the liquid phase properties at a temperature obtained as the average of boiling point and the ambient temperature). 


function [boilingTemperatureOne, boilingTemperatureTwo] = calBoilingPoint(liquidOne, liquidTwo, liquidMassFractionTwo, ambientPressure)

ambientPressure = ambientPressure / 101325;      % convert to bar

%% Choose the liquid 

water = strcmp(liquidOne, 'water');
rankWater = 1 ; 
if water ~= 1 
    water = strcmp(liquidTwo, 'water');
    rankWater = 2 ; 
end

urea = strcmp(liquidOne, 'urea');
rankUrea = 1 ; 
if urea ~= 1
    urea = strcmp(liquidTwo, 'urea');
    rankUrea = 2 ;
end

nHeptane = strcmp(liquidOne, 'nHeptane');
rankHeptane = 1 ; 
if nHeptane ~= 1
    nHeptane = strcmp(liquidTwo, 'nHeptane');
    rankHeptane = 2 ;
end

nDecane = strcmp(liquidOne, 'nDecane');
rankDecane = 1 ;
if nDecane ~= 1
    nDecane = strcmp(liquidTwo, 'nDecane');
    rankDecane = 2 ;
end

%% Calculate the boiling points

if water == 1           % Antoine Equation Parameters obtained from Properties of gases and liquids by Poling, Prausnitz and O'Connel
    
    A = 5.11564 ; 
    B = 1687.537 ;
    C = 230.17 ;
    
    boilingTemperature =  B / (A - log(ambientPressure)) - C + 273.15 ;
    
    if rankWater == 2 
        boilingTemperatureTwo = boilingTemperature ; 
    else
        boilingTemperatureOne = boilingTemperature ; 
    end
    
end
    
if urea == 1           % using the same equation as that of water, so it uses the boiling point of water as boiling point of urea. 
    
    A = 5.11564 ; 
    B = 1687.537 ;
    C = 230.17 ;
    
    boilingTemperature =  B / (A - log(ambientPressure)) - C + 273.15 ;
    
    if rankUrea == 2 
        boilingTemperatureTwo = boilingTemperature ; 
    else
        boilingTemperatureOne = boilingTemperature ; 
    end
    
end
    
if nHeptane == 1      % Ambrose and Ghiassee equation parameters obtained from Properties of gases and liquids by Poling, Prausnitz and O'Connel
    
    a = -7.77404 ; 
    b = 1.8564 ;
    c = -2.8298 ;
    d = -3.5070 ;
    Tc = 540.15 ;
    Pc = 27.35 ;
    
    reducedAmbientPressure = ambientPressure / Pc ;
    
    boilingTemperature = 300 ; 
        
    Tr = boilingTemperature / Tc ; 
    
    tau = 1 - Tr ; 
    
    reducedVaporPressure = exp(( a * tau + b * tau^1.5 + c * tau^2.5 + d * tau^5 )/Tr);
    
    while abs(reducedVaporPressure - reducedAmbientPressure)/reducedAmbientPressure > 0.01
        
        if reducedVaporPressure < reducedAmbientPressure
            
            boilingTemperature = boilingTemperature + 1 ; 
            
        else
            
            boilingTemperature = boilingTemperature - 1 ;
            
        end           
        
        Tr = boilingTemperature / Tc ; 
    
        tau = 1 - Tr ; 
    
        reducedVaporPressure = exp(( a * tau + b * tau^1.5 + c * tau^2.5 + d * tau^5 )/Tr);
        
    end
    
    if rankHeptane == 2 
        boilingTemperatureTwo = boilingTemperature ; 
    else
        boilingTemperatureOne = boilingTemperature ; 
    end
   
end

if nDecane == 1    % Ambrose and Ghiassee equation parameters obtained from Properties of gases and liquids by Poling, Prausnitz and O'Connel
   
    a = -8.60643 ; 
    b = 2.44659 ;
    c = -4.2925 ;
    d = -3.9080 ;
    Tc = 617.75 ;
    Pc = 21.05 ;
    
    reducedAmbientPressure = ambientPressure / Pc ;
    
    boilingTemperature = 300 ;
        
    Tr = boilingTemperature / Tc ; 
    
    tau = 1 - Tr ; 
    
    reducedVaporPressure = exp(( a * tau + b * tau^1.5 + c * tau^2.5 + d * tau^5 )/Tr);
    
    while abs(reducedVaporPressure - reducedAmbientPressure)/reducedAmbientPressure > 0.01
        
        if reducedVaporPressure < reducedAmbientPressure
            
            boilingTemperature = boilingTemperature + 1 ; 
            
        else
            
            boilingTemperature = boilingTemperature - 1 ;
            
        end           
        
        Tr = boilingTemperature / Tc ; 
    
        tau = 1 - Tr ; 
    
        reducedVaporPressure = exp(( a * tau + b * tau^1.5 + c * tau^2.5 + d * tau^5 )/Tr);
        
    end
    
    if rankDecane == 2 
        boilingTemperatureTwo = boilingTemperature ;
    else
        boilingTemperatureOne = boilingTemperature ; 
    end    
    
end

if le (liquidMassFractionTwo, 0)
        
    boilingTemperatureTwo = 0 ;
        
end


 

