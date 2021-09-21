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

%% This function calculates the reduced saturation pressure. 
%  Used for evaluating the density of a gas, using Peng-Robinson's EOS

function [reducedSatPressure] = calReducedSatPressure(nameOfGas, calculationMethod)

%% select the liquid
nitrogen = strcmp(nameOfGas,'N2');
air = strcmp(nameOfGas,'air');
waterVapor = strcmp(nameOfGas,'waterVapor');
ureaVapor = strcmp(nameOfGas,'ureaVapor');
nHeptaneVapor = strcmp(nameOfGas,'nHeptaneVapor');
nDecaneVapor = strcmp(nameOfGas,'nDecaneVapor');

%% select the correlation
mikhil = strcmp (calculationMethod, 'mikhil');
other = strcmp (calculationMethod, 'other');
abramzon = strcmp (calculationMethod, 'abramzon');
daif = strcmp (calculationMethod, 'daif');
sazhin = strcmp (calculationMethod, 'sazhin');
modifiedSazhin = strcmp(calculationMethod, 'modifiedSazhin');

%% Calculate the saturation vapor pressure based on selected correlation

if mikhil == 1 || abramzon == 1 || other == 1 || daif == 1 || sazhin == 1 || modifiedSazhin == 1
    
    if nitrogen == 1  || air == 1
        
        % Wagner Equation (7-3.3) for vapor pressure (Ambrose and Ghiassee, 1987) from the properties of gases and liquids by Poling, Prausnitz and O'Connel
        
        Tr = 0.7;
        tau = 1 - Tr ; 
    
        a = -6.1102 ;
        b = 1.2189 ;
        c = -0.69366 ;
        d = -1.89893 ;
        
        reducedSatPressure = exp((a * tau + b * tau^1.5 + c * tau^3 + d * tau^6 )/Tr) ;
  
    end

    if waterVapor == 1  
        
        % Wagner Equation (7-3.3) for vapor pressure (Ambrose and Ghiassee, 1987) from the properties of gases and liquids by Poling, Prausnitz and O'Connel
        
        Tr = 0.7;
        tau = 1 - Tr ; 
    
        a = -7.77224 ;
        b = 1.45684 ;
        c = -2.71942 ;
        d = -1.41336 ;
        
        reducedSatPressure = exp((a * tau + b * tau^1.5 + c * tau^3 + d * tau^6 )/Tr) ;
  
    end

    if nHeptaneVapor == 1
        
        % Wagner Equation (7-3.3) for vapor pressure (Ambrose and Ghiassee, 1987) from the properties of gases and liquids by Poling, Prausnitz and O'Connel
        
        Tr = 0.7 ;
        tau = 1 - Tr ; 
    
        a = -7.77404  ;
        b = 1.85614 ;
        c = -2.8298 ;
        d = -3.5070 ;
        
        reducedSatPressure = exp ((a * tau + b * tau^1.5 + c * tau^2.5 + d * tau^5 )/Tr);
    
    end
    
    if nDecaneVapor == 1
    
        % Wagner Equation (7-3.3) for vapor pressure (Ambrose and Ghiassee, 1987) from the properties of gases and liquids by Poling, Prausnitz and O'Connel
        
        Tr = 0.7 ;
        tau = 1 - Tr ; 
    
        a = -8.60643 ;
        b = 2.44659 ;
        c = -4.2925 ;
        d = -3.9080 ;
        
        reducedSatPressure = exp ((a * tau + b * tau^1.5 + c * tau^2.5 + d * tau^5 )/Tr); 

    end
    
end