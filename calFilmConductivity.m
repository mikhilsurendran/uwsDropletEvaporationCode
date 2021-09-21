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

%% This function calculates the combined conductivty of the constituent vapors of the film around the evaporating droplet.
%  Wilke's mixture rule is used to evaluate the conductivity of the film. Method of implementation can be improved using matrices. 

function [filmConductivity] = calFilmConductivity(K1, K2, K3, mu1, mu2, mu3, moleFraction1, moleFraction2, moleFraction3, molecularWeightOne, molecularWeightTwo, gas, calMethod)

%% select the correlation
mikhil = strcmp(calMethod, 'mikhil');

%% select the ambient gas
nitrogen = strcmp(gas,'N2');
air = strcmp(gas,'air');

if air == 1
    
    gasMolecularWeight = 28.01 ; % g/mol
    
elseif nitrogen == 1
    
    gasMolecularWeight = 28.97 ; % g/mol
    
end

%% Calculation of the conductivity of the film

if mikhil==1           % Using Wilke Mixing Rule
    
    filmConductivity = 0 ;
    
    moleFrac = zeros(1,3);
    moleFrac(1,1) = moleFraction1 ; 
    moleFrac(1,2) = moleFraction2 ;
    moleFrac(1,3) = moleFraction3 ; 

    mu = zeros(1,3);
    mu(1,1) = mu1 ; 
    mu(1,2) = mu2 ;
    mu(1,3) = mu3 ;
    
    K = zeros (1,3);
    K(1,1) = K1 ; 
    K(1,2) = K2 ;
    K(1,3) = K3 ;
    
    % for bi component liquids
    if mu(1,2) ~= 0
        
        phiWilke = zeros (3,3) ; 
        phiWilke(1,1) = 1 ; 
        phiWilke(1,2) = ( ( 1 + (mu(1,1)/mu(1,2))^(0.5) * (molecularWeightTwo/molecularWeightOne)^(0.25) ) ^ 2 ) / (( 8^(0.5) ) * (1 + molecularWeightOne/molecularWeightTwo)^(0.5)) ;  
        phiWilke(1,3) = ( ( 1 + (mu(1,1)/mu(1,3))^(0.5) * (gasMolecularWeight/molecularWeightOne)^(0.25) ) ^ 2 ) / (( 8^(0.5) ) * (1 + molecularWeightOne/gasMolecularWeight)^(0.5)) ; 
        phiWilke(2,1) = ( ( 1 + (mu(1,2)/mu(1,1))^(0.5) * (molecularWeightOne/molecularWeightTwo)^(0.25) ) ^ 2 ) / (( 8^(0.5) ) * (1 + molecularWeightTwo/molecularWeightOne)^(0.5)) ; 
        phiWilke(2,2) = 1 ; 
        phiWilke(2,3) = ( ( 1 + (mu(1,2)/mu(1,3))^(0.5) * (gasMolecularWeight/molecularWeightTwo)^(0.25) ) ^ 2 ) / (( 8^(0.5) ) * (1 + molecularWeightTwo/gasMolecularWeight)^(0.5)) ; 
        phiWilke(3,1) = ( ( 1 + (mu(1,3)/mu(1,1))^(0.5) * (molecularWeightOne/gasMolecularWeight)^(0.25) ) ^ 2 ) / (( 8^(0.5) ) * (1 + gasMolecularWeight/molecularWeightOne)^(0.5)) ; 
        phiWilke(3,2) = ( ( 1 + (mu(1,3)/mu(1,2))^(0.5) * (molecularWeightTwo/gasMolecularWeight)^(0.25) ) ^ 2 ) / (( 8^(0.5) ) * (1 + gasMolecularWeight/molecularWeightTwo)^(0.5)) ; 
        phiWilke(3,3) = 1 ; 
        
        for i = 1:3
            
            numerator = (moleFrac(1,i) * K(1,i));
            
            denominator = moleFrac(1,i) ; 
            
            for j = 1:3 
                
                if i ~= j 
                    
                    denominator = denominator + ( moleFrac(1,j) * phiWilke(i,j) ) ; 
                    
                end
                
            end
            
            filmConductivity = filmConductivity + (numerator / denominator);
            
        end
    
    % for single component liquids
    else
        
        moleFrac = zeros(1,2);
        moleFrac(1,1) = moleFraction1 ;
        moleFrac(1,2) = moleFraction3 ; 
        
        mu = zeros(1,2);
        mu(1,1) = mu1 ; 
        mu(1,2) = mu3 ;
        
        K = zeros (1,2);
        K(1,1) = K1 ; 
        K(1,2) = K3 ;
        
        phiWilke = zeros (2,2) ;
        phiWilke(1,1) = 1 ; 
        phiWilke(1,2) = ( ( 1 + (mu(1,1)/mu(1,2))^(0.5) * (gasMolecularWeight/molecularWeightOne)^(0.25) ) ^ 2 ) / (( 8^(0.5) ) * (1 + molecularWeightOne/gasMolecularWeight)^(0.5)) ;  
        phiWilke(2,1) = ( ( 1 + (mu(1,2)/mu(1,1))^(0.5) * (molecularWeightOne/gasMolecularWeight)^(0.25) ) ^ 2 ) / (( 8^(0.5) ) * (1 + gasMolecularWeight/molecularWeightOne)^(0.5)) ; 
        phiWilke(2,2) = 1 ;
        
        for i = 1:2
            
            numerator = (moleFrac(1,i) * K(1,i));
            
            denominator = moleFrac(1,i) ; 
            
            for j = 1:2 
                
                if i ~= j 
                    
                    denominator = denominator + ( moleFrac(1,j) * phiWilke(i,j) ) ; 
                    
                end
                
            end
            
            filmConductivity = filmConductivity + (numerator / denominator);
            
        end
        
    end
    
end       
