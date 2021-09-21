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

%% This function calculates the partial pressure for the vapor mixture. 
%  By default droplet surface temperature is used to evaluate the
%  saturation pressure. 
%  Inputs and outputs are in SI units

function [satPressureOne, satPressureTwo, totalVaporPressure] = calPartialPressure(saturationPressureOne, saturationPressureTwo, ambientTemperature, surfaceTemperature, liqCompTwoMassFraction, componentOneBoilPt, filmRulePara, componentOneMolWt, componentTwoMolWt, satPressLaw, calMethod, isItUWS)

%% select the correlation
mikhil = strcmp (calMethod, 'mikhil');

%% select the partial pressure law 
raoults = strcmp (satPressLaw, 'raoults'); 
nrtl = strcmp (satPressLaw, 'NRTL');

%% Calculate the partial vapor pressure based on selected correlation

[moFracOne, moFracTwo] = calLiquidMoleFrac(liqCompTwoMassFraction, componentOneMolWt, componentTwoMolWt) ;

if mikhil == 1 
    
    if raoults == 1 
            
        satPressureOne = moFracOne * saturationPressureOne ;   
        satPressureTwo = moFracTwo * saturationPressureTwo ; 
        
        if strcmp(isItUWS,'urea')==1  && le(surfaceTemperature, 406.15)
            
            satPressureTwo = 0 ;    % urea does not evaporate / melt below 406.15 K 
            
        end 
        
        totalVaporPressure = satPressureOne + satPressureTwo ;
   
    end
    
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if nrtl == 1 
        
        if strcmp(isItUWS,'urea')==1 
            
            temperature = surfaceTemperature ; 
        
            c = [0; 0] ; 
            alpha = 0.3 ; 
            a = [7.659; -1.536] ;  
            b = [-1463; -56.67] ;
        
            tau = zeros(2,1);
        
            tau(1,1) = a(1,1) + b(1,1)/temperature + c(1,1)/temperature^2 ;
            tau(2,1) = a(2,1) + b(2,1)/temperature + c(2,1)/temperature^2 ;
        
            G = zeros(2,1);
        
            G(1,1) = exp(-alpha * tau(1,1)) ;
            G(2,1) = exp(-alpha * tau(2,1)) ;
        
            gammaW = exp( moFracTwo^2 * ( tau(2,1) * ( G(2,1)/(moFracOne + moFracTwo * G(2,1)) )^2 + (tau(1,1)*G(1,1)/(moFracTwo + moFracOne * G(1,1))^2) ) ) ;
            gammaU = exp( moFracOne^2 * ( tau(1,1) * ( G(1,1)/(moFracTwo + moFracOne * G(1,1)) )^2 + (tau(2,1)*G(2,1)/(moFracOne + moFracTwo * G(2,1))^2) ) ) ;
    
            satPressureOne = moFracOne * gammaW * saturationPressureOne ;
            satPressureTwo = moFracTwo * gammaU * saturationPressureTwo ;
        
            if le(surfaceTemperature, 406.15)
            
                satPressureTwo = 0 ;    % urea does not evaporate / melt below 406.15 K 
            
            end
        
            totalVaporPressure = satPressureOne + satPressureTwo ;
            
        else
            
            msg = 'In this code NRTL has been implemented for UWS alone. Please choose another partial pressure law. '; 
            error(msg) ; 
            
        end
    
    end
    
end