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

%% This function calculates the mass of the bi-component droplet. 

function [massOfDrop, compOneMass, compTwoMass, dMassDt, dMassOneDt, dMassTwoDt] = calDropletMass(dropDia, dropTemp, compOneBP, compTwoBP, compTwoMF, liquidOne, liquidTwo, ambPress, ambTemp, dTemperatureDt, dDiameterDt, dMassFractionDt, calculationMethod, character, exptOrCorrln)

FRP = 0 ; % film rule parameter

dropVolume = pi() * dropDia .^3 / 6 ;

[numRows, numCols] = size(dropDia) ;

massOfDrop = zeros(numRows,numCols);
compTwoMass = zeros(numRows,numCols);
dVolDt = zeros(numRows,numCols);
dMassDt = zeros(numRows,numCols);
dMassOneDt = zeros(numRows,numCols);
dMassTwoDt = zeros(numRows,numCols); 

for count = 1:numRows
    
    %% Evaluate the density and mass at all the time steps
    [densityOne, dOnedT] = calRhoLiquid(ambPress, ambTemp, dropTemp(1,1), compOneBP, dropTemp(count,1), FRP, liquidOne, liquidTwo, compTwoMF(count,1), calculationMethod, character, exptOrCorrln);
    [densityTwo, dTwodT] = calRhoLiquid(ambPress, ambTemp, dropTemp(1,1), compTwoBP, dropTemp(count,1), FRP, liquidTwo, liquidTwo, compTwoMF(count,1), calculationMethod, character, exptOrCorrln);
    
    densityMix = calMixDensity(densityOne, densityTwo, compTwoMF(count,1), calculationMethod) ;
    
    massOfDrop(count,1) = densityMix * dropVolume(count,1);
    compTwoMass(count,1) = massOfDrop(count,1) * compTwoMF(count,1);
    
    %% Evaluate the rate of change of liquid density at all time steps 
    % Mixture density is a function of temperature and liquid mass fraction
    % Both temperature and mass fraction are time varying.
    % Hence the rate of change of liquid density at all time steps is evaluated using the derivative of density with respect to temperature (dOnedT, dTwodT), and the rate of change of mass fractions
    
    dOneBydTime = dOnedT * dTemperatureDt(count,1) ; 
    dTwoBydTime = dTwodT * dTemperatureDt(count,1) ; 
    
    % calculate dRhodt
    firstTerm = ( (1-compTwoMF(count,1))/densityOne^2 ) * dOneBydTime ;
    if compTwoMF > 0 
        secondTerm = ( compTwoMF(count,1)/densityTwo^2 ) * dTwoBydTime ;
        thirdTerm = ( 1/densityOne - 1/densityTwo ) * (dMassFractionDt(count,1)) ;
    else
        secondTerm = 0 ; 
        thirdTerm = 0 ; 
    end

    dDensityDt = (densityMix)^2 * (firstTerm + secondTerm + thirdTerm ) ;
    
    %% Evaluate the rate of change of droplet Volume at all time steps
    dVolDt(count,1) = (pi/6) * ( 3 * dropDia(count,1)^2 ) * (dDiameterDt(count,1)) ;
     
    %% Evaluate the rate of change of mass at all time steps
    dMassDt(count,1) = dropVolume(count,1) * dDensityDt + densityMix * dVolDt(count,1) ;    % m = rho*V, dmdt = rho * dVdt + V * dRhodt 
    dMassTwoDt(count,1) = compTwoMF(count,1) * dMassDt(count,1) + massOfDrop(count,1) * dMassFractionDt(count,1) ; % m2 = m*Y2, dm2dt = m * dY2dt + Y2 * dmdt
    dMassOneDt(count,1) = dMassDt(count,1) - dMassTwoDt(count,1) ; % dm1dt = dmdt - dm2dt 
    
end

%% Convert Mass from kg to micro grams
massOfDrop = massOfDrop * 1e9;    % converting from kg to ug   
compTwoMass = compTwoMass * 1e9;  % converting from kg to ug
compOneMass = massOfDrop - compTwoMass;  % converting from kg to ug 

dMassDt = -1 * dMassDt * 1e9 ; % converting from kg to ug (multiplied by -1 because it is mass loss)
dMassOneDt = -1 * dMassOneDt * 1e9 ; % converting from kg to ug (multiplied by -1 because it is mass loss)
dMassTwoDt = -1 * dMassTwoDt * 1e9 ; % converting from kg to ug (multiplied by -1 because it is mass loss)

