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

%% This function calculates the viscosity of the gas/vapour. It may be suitably modified to use other expressions. 

function [viscosity] = calViscosity(ambientPressure, ambientTemperature, surfaceTemperature, filmRulePara, gas, calMethod)

%% select the correlation
mikhil = strcmp(calMethod, 'mikhil');

%% select the fluid
nitrogen = strcmp(gas,'N2');
air = strcmp(gas,'air');
oxygen = strcmp(gas,'O2');
argon = strcmp(gas,'Ar');
waterVapor = strcmp(gas,'waterVapor');
ureaVapor = strcmp(gas,'ureaVapor');
nHeptaneVapor = strcmp(gas,'nHeptaneVapor');
nDecaneVapor = strcmp(gas,'nDecaneVapor');

%% calculate the dynamic viscosity

if mikhil == 1
    
    filmTemperature = surfaceTemperature*(1-filmRulePara)+ ambientTemperature*filmRulePara;
    
    if nDecaneVapor == 1              % using correlation from abramzon's paper
        
        viscosity = (0.564 + 1.75 * 10 ^(-3) * (filmTemperature - 300)) * 10 ^(-5) ; 
        
    elseif nHeptaneVapor == 1         % using correlation from abramzon's paper
        
        viscosity = (0.620 + 1.80 * 10 ^(-3) * (filmTemperature - 300)) * 10 ^(-5) ; 
        
    else                              % not from the original paper, but using correlations from Chemical Properties Handbook by Carl L Yaws
        
        if nitrogen == 1
    
            A = 42.606;
            B = 4.75*10^(-1);
            C = -9.88*10^(-5);
    
        end

        if air == 1                   % not configured. 
    
            A = 0;
            B = 0;
            C = 0;
    
        end
        
        if argon == 1
    
            A = 44.997;
            B = 6.3892*10^(-1);
            C = -1.2455*10^(-4);
    
        end
    
        if oxygen == 1
    
            A = 44.224;
            B = 5.62*10^(-1);
            C = -1.13*10^(-4);
    
        end

        if waterVapor == 1
    
            A = -36.826;
            B = 4.29*10^(-1);
            C = -1.62*10^(-5);
    
        end

        if ureaVapor == 1
   
            if ge(surfaceTemperature,406.15)
           
                A = -13.895;
                B = 2.7802*10^(-1);
                C = -3.8420*10^(-5);
        
            else
        
                A = 0;
                B = 0;
                C = 0;
        
            end
    
        end

        viscosity = ( A + B * filmTemperature + C * filmTemperature^2 )*10^(-7);
    
    end
    
end
