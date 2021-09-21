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

%% This function evaluates the molefractions of the constituents of the ambient based on their mass fractions and molecular weights.

function [compOneMoleFraction, compTwoMoleFraction, ambGasMoleFraction] = calVapMoleFrac_ambient(compOneMassFraction, compTwoMassFraction, molecularWeightOne, molecularWeightTwo, ambMolecularWeight)

ambGasMassFraction = 1 - (compOneMassFraction + compTwoMassFraction);

if compOneMassFraction > 0 && compTwoMassFraction > 0
    
    compOneMoleFraction = 1/( 1 + (compTwoMassFraction/molecularWeightTwo) * (molecularWeightOne/compOneMassFraction) + (ambGasMassFraction/ambMolecularWeight) * (molecularWeightOne/compOneMassFraction) );

    compTwoMoleFraction = 1/( 1 + (compOneMassFraction/molecularWeightOne) * (molecularWeightTwo/compTwoMassFraction) + (ambGasMassFraction/ambMolecularWeight) * (molecularWeightTwo/compTwoMassFraction) );

    ambGasMoleFraction = 1 - (compOneMoleFraction + compTwoMoleFraction); 
    
end

if compOneMassFraction == 0 && compTwoMassFraction > 0
    
    compOneMoleFraction = 0 ;

    compTwoMoleFraction = 1/( 1 + (ambGasMassFraction/ambMolecularWeight) * (molecularWeightTwo/compTwoMassFraction) );

    ambGasMoleFraction = 1 - (compOneMoleFraction + compTwoMoleFraction); 
    
end

if compTwoMassFraction == 0 && compOneMassFraction > 0
    
    compTwoMoleFraction = 0 ;

    compOneMoleFraction = 1/( 1 + (ambGasMassFraction/ambMolecularWeight) * (molecularWeightOne/compOneMassFraction) );

    ambGasMoleFraction = 1 - (compOneMoleFraction + compTwoMoleFraction); 
    
end

if compTwoMassFraction == 0 && compOneMassFraction == 0
    
    compTwoMoleFraction = 0 ;

    compOneMoleFraction = 0 ;

    ambGasMoleFraction = 1 ; 
    
end

if compOneMoleFraction < 1e-5
    compOneMoleFraction = 0 ; 
elseif compOneMoleFraction > 0.99999
    compOneMoleFraction = 1 ; 
end

if compTwoMoleFraction < 1e-5
    compTwoMoleFraction = 0 ; 
elseif compTwoMoleFraction > 0.99999
    compTwoMoleFraction = 1 ; 
end

if ambGasMoleFraction < 1e-5
    ambGasMoleFraction = 0 ; 
elseif ambGasMoleFraction > 0.99999
    ambGasMoleFraction = 1 ; 
end
