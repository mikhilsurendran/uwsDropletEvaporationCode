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

%% This function calculates the thermal conductivity of the gas/vapour. It may be suitably modified to use other expressions. 
%  Another function calVaporFilmConductivity(...) is used to calculate the effective conductivity of the
%  film surrounding the droplet using individual vapor conductivities calculated
%  using calVaporConductivity(...)

function [conductivity] = calVaporConductivity(ambientPressure, ambientTemperature, surfaceTemperature, filmRulePara, gas, calMethod)

%% select the correlation
mikhil = strcmp(calMethod, 'mikhil');

%% select the fluid
nitrogen = strcmp(gas,'N2');
air = strcmp(gas,'air');
argon = strcmp(gas,'Ar');
oxygen = strcmp(gas,'O2');
waterVapor = strcmp(gas,'waterVapor');
ureaVapor = strcmp(gas,'ureaVapor');
nHeptaneVapor = strcmp(gas,'nHeptaneVapor');
nDecaneVapor = strcmp(gas,'nDecaneVapor');

%% Calculate vapor conductivity

if mikhil == 1
    
    filmTemperature = surfaceTemperature*(1-filmRulePara)+ ambientTemperature*filmRulePara;
    
    
    if nDecaneVapor == 1         % using correlations from abramzon's paper
       
        conductivity = 2.9 * 10 ^(-6) * (filmTemperature/300) ^ (1.8) ;    % kcal / m s K 
        
        conductivity = conductivity * 4.184 * 1000 ; % W / m K
        
        
    elseif nHeptaneVapor == 1    % using correlations from abramzon's paper 
    
        conductivity = 2.9 * 10 ^(-6) * (filmTemperature/300) ^ (1.8) ;    % kcal / m s K 
        
        conductivity = conductivity * 4.184 * 1000 ; % W / m K
        
         
    else                         % using correlations given in Chemical Properties Handbook by Carl L Yaws.

        if nitrogen == 1
    
            A = 0.00309;
            B = 7.5930*10^-5;
            C = -1.1014*10^-8;
    
        end

        if air == 1    % not configured
    
            A = 0;
            B = 0;
            C = 0;
    
        end
        
        if argon == 1
    
            A = 0.00548;
            B = 4.3869*10^(-5);
            C = -6.8141*10^(-9);
    
        end
    
        if oxygen == 1
    
            A = 0.00121;
            B = 8.6157*10^(-5);
            C = -1.3346*10^(-9);
    
        end

        if waterVapor == 1
    
            A = 0.00053;
            B = 4.7093*10^-5;
            C = 4.9551*10^-8;
    
        end

        if ureaVapor == 1
    
            if ge(surfaceTemperature, 406.15)  % urea vaporization occurs at temperatures above 133 degC (406.15 K)
        
                A = -0.01444;
                B = 7.2675 * 10^-5;
                C = 1.2011 * 10^-8;
        
            else
        
                A = 0;
                B = 0;
                C = 0; 
        
            end
    
        end
  
        conductivity = A + B*filmTemperature + C*filmTemperature^2;   
        
    end
    
end