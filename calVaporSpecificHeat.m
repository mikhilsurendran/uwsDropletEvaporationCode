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

%% This function calculates the specific heat of the gas/vapour. It may be suitably modified to use other expressions. 
%  Another function calVaporFilmSpecificHeat(...) is used to calculate the effective specific heat of the
%  film surrounding the droplet using individual specific heats calculated
%  using calVaporSpecificHeat(...)

function [specificHeat] = calVaporSpecificHeat(ambientPressure, ambientTemperature, surfaceTemperature, filmRulePara, gas, calMethod)

%% select the fluid
nitrogen = strcmp(gas,'N2');
air = strcmp(gas,'air');
waterVapor = strcmp(gas,'waterVapor');
ureaVapor = strcmp(gas,'ureaVapor');
nHeptaneVapor = strcmp(gas,'nHeptaneVapor');
nDecaneVapor = strcmp(gas,'nDecaneVapor');

%% select the correlation
mikhil = strcmp(calMethod, 'mikhil');

%% calculate vapor specific heats

if mikhil == 1 
    
    filmTemperature = surfaceTemperature*(1-filmRulePara) + ambientTemperature*filmRulePara;
    
    if nDecaneVapor == 1   

        A = 1.6720 * 10^5 ;  % from Perry's Handbook
        B = 5.3530 * 10^5 ;
        C = 1.6141 * 10^3 ;
        D = 3.7820 * 10^5 ;
        E = 742 ;
        
        molecularWeight = 142.285;
    
        specificHeat = A + B * ((C/filmTemperature)/sinh(C/filmTemperature))^2 + D * ((E/filmTemperature)/cosh(E/filmTemperature))^2 ; % J/kmol K
    
        specificHeat = specificHeat / (molecularWeight) ; % J / kg K
        
        
    else                           % not using correlations from the original paper, but from those given in Perrys Handbook
        
        if nitrogen == 1
    
            A = 0.2911 * 10^5 ;  % from Perry's Handbook
            B = 0.0861 * 10^5 ;
            C = 1.7016 * 10^3 ;
            D = 0.0010 * 10^5 ;
            E = 909.79 ;
    
            molecularWeight = 28.01; % kg/kmol
    
            specificHeat = A + B * ((C/filmTemperature)/sinh(C/filmTemperature))^2 + D * ((E/filmTemperature)/cosh(E/filmTemperature))^2 ; % J/kmol K
    
            specificHeat = specificHeat / (molecularWeight) ; % J / kg K
    
        elseif air == 1
    
            A = 0.2896 * 10^5 ;  % from Perry's Handbook
            B = 0.0939 * 10^5 ;
            C = 3.0120 * 10^3 ;
            D = 0.0758 * 10^5 ;
            E = 1484 ;
    
            molecularWeight = 28.97; % g/mol
        
            specificHeat = A + B * ((C/filmTemperature)/sinh(C/filmTemperature))^2 + D * ((E/filmTemperature)/cosh(E/filmTemperature))^2 ; % J/kmol K
    
            specificHeat = specificHeat / (molecularWeight) ; % J / kg K
    

        elseif waterVapor == 1
    
            A = 0.3336 * 10^5 ;  % from Perry's Handbook
            B = 0.2679 * 10^5 ;
            C = 2.6105 * 10^3 ;
            D = 0.0890 * 10^5 ;
            E = 1169 ;
    
            molecularWeight = 18.02; % g/mol
    
            specificHeat = A + B * ((C/filmTemperature)/sinh(C/filmTemperature))^2 + D * ((E/filmTemperature)/cosh(E/filmTemperature))^2 ; % J/kmol K
    
            specificHeat = specificHeat / (molecularWeight) ; % J / kg K
    
        elseif ureaVapor == 1
    
            if ge(surfaceTemperature, 406.15) % urea vaporization occurs at temperatures above 133 degC (406.15 K)
        
                A = 24.856;
                B = 1.4437*10^(-1);
                C = 3.8088*10^(-5);
                D = -1.1007*10^(-7);
                E = 3.9161*10^(-11);
        
                molecularWeight = 60.06; % g/mol
        
                specificHeat = A + B*filmTemperature + C*filmTemperature^2 + D*filmTemperature^3 + E*filmTemperature^4;   %J/mol K
    
                specificHeat = specificHeat / (molecularWeight / 1000) ; % J / kg K
        
            else
    
                specificHeat = 0 ; % J / kg K
        
            end
    
        elseif nHeptaneVapor == 1
    
            A = 1.2015 * 10^5 ;  % from Perry's Handbook
            B = 4.0010 * 10^5 ;
            C = 1.6766 * 10^3 ;
            D = 2.7400 * 10^5 ;
            E = 756.4 ;
        
            molecularWeight = 100.204;
    
            specificHeat = A + B * ((C/filmTemperature)/sinh(C/filmTemperature))^2 + D * ((E/filmTemperature)/cosh(E/filmTemperature))^2 ; % J/kmol K
    
            specificHeat = specificHeat / (molecularWeight) ; % J / kg K       
        
        end
        
    end
    
end