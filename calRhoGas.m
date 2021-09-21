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

%% This function calculates the density of the gas. It may be suitably modified to use other laws. 

function [density] = calRhoGas(ambientPressure, ambientTemperature, surfaceTemperature, filmRulePara, gas, rule, calMethod)

%% select the correlation
mikhil = strcmp(calMethod, 'mikhil');

%% select the fluid
nitrogen = strcmp(gas,'N2');
air = strcmp(gas,'air');
waterVapor = strcmp(gas,'waterVapor');
ureaVapor = strcmp(gas,'ureaVapor');
nHeptaneVapor = strcmp(gas,'nHeptaneVapor');
nDecaneVapor = strcmp(gas,'nDecaneVapor');

%% calculate the density

if mikhil == 1

    filmTemperature = surfaceTemperature*(1-filmRulePara) + ambientTemperature*filmRulePara;

    if nitrogen == 1
    
        molecularWeight = 28.01 ; % g/mol
        Tc = getCriticalTemperature(gas);
        Pc = getCriticalPressure(gas);
    
    elseif air == 1
    
        molecularWeight = 28.97 ; % g/mol
        Tc = getCriticalTemperature(gas);
        Pc = getCriticalPressure(gas);

    elseif waterVapor == 1
    
        molecularWeight = 18.02 ; % g/mol
        Tc = getCriticalTemperature(gas);
        Pc = getCriticalPressure(gas);
    
    elseif ureaVapor == 1
        
        molecularWeight = 60.06 ; % g/mol 
    
    elseif nHeptaneVapor == 1
    
        molecularWeight = 100.204 ; % g/mol
        Tc = getCriticalTemperature(gas);
        Pc = getCriticalPressure(gas);
    
    elseif nDecaneVapor == 1
    
        molecularWeight = 142.285 ; % g/mol
        Tc = getCriticalTemperature(gas);
        Pc = getCriticalPressure(gas);
    
    end   
    
    % density is calculated based on the selected gas law
    ideal = strcmp(rule,'ideal');
    real = strcmp(rule,'real');

    if ideal == 1
    
        R = 8.314 ; % J/K mol
    
        density = ambientPressure*molecularWeight/(1000*R*filmTemperature); % kg/m^3
    
        if ureaVapor == 1 && surfaceTemperature < 406.15  % urea vaporization occurs at temperatures above 133 degC (406.15 K) 
        
            density = 0 ;
        
        end
    
    end

    if real == 1   % Peng - Robinson's equation of state
        
        filmTemperature = surfaceTemperature * (1-filmRulePara) + ambientTemperature * filmRulePara ;
        
        if ureaVapor == 1  % density of urea vapor is evaluated using ideal gas equation 
            
            density = ambientPressure*molecularWeight/(1000*R*filmTemperature); % kg/m^3
            
            if surfaceTemperature < 406.15  % urea vaporization occurs at temperatures above 133 degC (406.15 K) 
        
                density = 0 ;
                
            end
            
        else
            
            T = filmTemperature ; 
        
            R = 8.314 ; % J/K mol
        
            PrSat = calReducedSatPressure (gas, calMethod) ; 
        
            omega = - log10(PrSat) - 1 ;        % Acentric factor
        
            Tr = T / Tc  ; 
        
            k = 0.37464 + 1.5422 * omega - 0.26922 * omega^2 ; 
        
            a = 0.45724 * (R * Tc)^2 / Pc * ( 1 + k * ( 1 - Tr^(0.5) ) )^2 ; 
        
            b = 0.07780 * (R * Tc)/Pc ; 
        
            molarVolume = ( R * T ) / ambientPressure ;      % m^3 / mole (using ideal gas law)  
        
            Vm = molarVolume / 10 ;                      % initialGuess 
        
            increment = Vm / 100 ;
        
            pressure = 0 ; 
        
            while abs ((pressure - ambientPressure)/pressure) > 0.001 
            
                pressure = ( R * T ) / (Vm - b) -  a / (Vm^2 + 2*Vm*b - b^2) ;
            
                Vm = Vm + increment ; 
            
            end
        
            density = molecularWeight / ( Vm * 1000 ) ; % kg / m^3 
        
        end

    end
    
end