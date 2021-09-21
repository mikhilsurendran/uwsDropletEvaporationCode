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

%% This function evaluates and returns various thermophysical properties and non-dimensional numbers.
%  This fucntion is called only after the solution is complete, and the output
%  from this function is written to text files. 
%  This function is identical to evaluateProperties.m, but the number of
%  parameters returned to the calling function are higher here. 
%  While evaluateProperties.m returns only those properties and non-dimensional numbers needed to define
%  the derivatives of the variables being solved for, this functions
%  returns all the thermophysical properties and non-dimensional numbers, which are then saved as output in text files.  

function [compOneLiquidMoFr, compTwoLiquidMoFr, rhoLiqOne, dRhoBydTempOne, rhoLiqTwo, dRhoBydTempTwo, rhoLiq, cpLiqOne, cpLiqTwo, cpLiq, satPressCompOne, satPressCompTwo, totalSatPressure, compOneVapMoF_s, compTwoVapMoF_s, vapMoF_s, ambGasMoF_s, compOneVapMF_s, compTwoVapMF_s, ambGasMF_s, compOneVapMoF_inf, compTwoVapMoF_inf, ambGasMoF_inf, compOneVapMF_film, compTwoVapMF_film, ambGasMF_film, compOneVapMoF_film, compTwoVapMoF_film, ambGasMoF_film, filmMixMW, cpAmbGas, cpCompOneVap, cpCompTwoVap, cpFilm, condOneVap, condTwoVap, condN2, condO2, condAr, condAir, condAmbGas, condVapFilm, binDiffusivityOne, binDiffusivityTwo, binDiffusivityTwoInOne, muN2FilmTemp, muO2FilmTemp, muArFilmTemp, muAirFilmTemp, muAmbGas, muVapOne, muVapTwo, muN2AmbTemp, muVapFilm, rhoCompOneVap, rhoCompTwoVap, rhoInf, rhoAir, rhoVapFilm, compOneEnth, compTwoEnth, mixEnthalpy, cpVapor, LeCompOne, LeCompTwo, Re, Pr, ScCompOne, ScCompTwo, Sh0CompOne, Sh0CompTwo, BmCompOne, BmCompTwo, ShStarCompOne, ShStarCompTwo, Nu0, Bt, phi, NuStaar] = evaluateProperties_postSolvingTheODE(dropSurfTemp, dropDia, compTwoLiqMaFr, relativeVel)

intermediateSaveOption = 'no' ; 

% Copy paste the contents from the function evaluateProperties.m below this
% line. The functions used to evaluate the properties are exactly the same,
% but the number of parameters returned are different. 

%% Read user inputs from file

inputFile ;           % read the user inputs from file

iterationNumber ;     % read the iteration number from file

boilingPoints ;       % read the boiling points from file

vaporOne = strcat(componentOne.Name,'Vapor');
vaporTwo = strcat(componentTwo.Name,'Vapor');
vaporMix = 'vaporMix';

mikhil = strcmp(method, 'mikhil');                      % select the method specified by the user.

naturalConvection = strcmp(isNaturalConvectionIncluded, 'yes');     % check whether to include natural convection effects or not in case the relative air velocity is zero.

%% Calculate the liquid mole fractions

[compOneLiquidMoFr, compTwoLiquidMoFr] = calLiquidMoleFrac(compTwoLiqMaFr, componentOne.MW, componentTwo.MW);

%% Calculate the other liquid phase properties

[rhoLiqOne, dRhoBydTempOne, rhoLiqTwo, dRhoBydTempTwo, rhoLiq, cpLiqOne, cpLiqTwo, cpLiq, satPressCompOne, satPressCompTwo, totalSatPressure] = calLiquidPhaseProperties(ambient.Pressure, ambient.Temperature, droplet.InitialSurfaceTemp, dropSurfTemp, compTwoLiqMaFr, compTwoLiquidMoFr, compOneLiquidMoFr, componentOne.MW, componentTwo.MW, iterNo, intermediateSaveOption, componentOne.Name, componentTwo.Name, compOneBoilingPoint, compTwoBoilingPoint, partialPressureLaw, method, liquidPropertyCharacter, uwsDensityCalculationMethod, ureaVaporPressureCorrelation, componentTwoDepletionCalculationMethod); 

%% Calculate vapor phase properties

[compOneVapMoF_s, compTwoVapMoF_s, vapMoF_s, ambGasMoF_s, compOneVapMF_s, compTwoVapMF_s, ambGasMF_s, compOneVapMoF_inf, compTwoVapMoF_inf, ambGasMoF_inf, compOneVapMF_film, compTwoVapMF_film, ambGasMF_film, compOneVapMoF_film, compTwoVapMoF_film, ambGasMoF_film, filmMixMW, cpAmbGas, cpCompOneVap, cpCompTwoVap, cpFilm, condOneVap, condTwoVap, condN2, condO2, condAr, condAir, condAmbGas, condVapFilm, binDiffusivityOne, binDiffusivityTwo, binDiffusivityTwoInOne, muN2FilmTemp, muO2FilmTemp, muArFilmTemp, muAirFilmTemp, muAmbGas, muVapOne, muVapTwo, muN2AmbTemp, muVapFilm, rhoCompOneVap, rhoCompTwoVap, rhoInf, rhoAir, rhoVapFilm, compOneEnth, compTwoEnth, mixEnthalpy] = calVaporPhaseProperties(ambient.Pressure, ambient.Temperature, dropSurfTemp, satPressCompOne, satPressCompTwo, compTwoLiqMaFr, ambient.ComponentOneMF, ambient.ComponentTwoMF, componentOne.MW, componentTwo.MW, ambient.Gas, ambient.GasMW, vaporOne, vaporTwo, equationOfState, iterNo, intermediateSaveOption, method, diffusivityModel, isDissociationEnergyIncluded, ureaVaporizationEnthalpyCorrelation, componentTwoDepletionCalculationMethod);

if compTwoVapMoF_s > 0
    
    cpVapor = cpCompOneVap * compOneVapMF_s / (compOneVapMF_s + compTwoVapMF_s)  + cpCompTwoVap * compTwoVapMF_s / (compOneVapMF_s + compTwoVapMF_s) ; 
    
else
    
    cpVapor = cpCompOneVap ; 
    
end

%% Lewis number

LeCompOne = calLewis(condVapFilm, condAmbGas, cpFilm, cpAmbGas, rhoCompOneVap, rhoCompTwoVap, rhoAir, rhoVapFilm, binDiffusivityOne, ambient.Temperature, dropSurfTemp, iterNo, vaporOne, intermediateSaveOption, method); 

if compTwoVapMoF_s > 0
    LeCompTwo = calLewis(condVapFilm, condAmbGas, cpFilm, cpAmbGas, rhoCompOneVap, rhoCompTwoVap, rhoAir, rhoVapFilm, binDiffusivityTwo, ambient.Temperature, dropSurfTemp, iterNo, vaporTwo, intermediateSaveOption, method);
else
    LeCompTwo = 0 ; 
end

%% Reynolds number 

Re = calReynolds(rhoInf, relativeVel, dropDia, muVapFilm, ambient.Temperature, iterNo, intermediateSaveOption, method);

%% Grashof number 

if naturalConvection == 1 
    
    Gr = calGrashof(rhoInf,  muVapFilm, dropDia, ambient.Temperature, dropSurfTemp, iterNo, intermediateSaveOption, method);
    
else 
    
    Gr = 0 ; 
    
end

%% Prandtl number 

Pr = calPrandtl(muVapFilm, cpFilm, condVapFilm, ambient.Temperature, iterNo, intermediateSaveOption, method);

%% Nusselt number 

Nu0 = calNusselt(Re, Gr, Pr, ambient.Temperature, iterNo, intermediateSaveOption, nusseltCorrelation);

%% Schmidt number 

ScCompOne = calSchmidt(muVapFilm, rhoVapFilm, binDiffusivityOne, ambient.Temperature, dropSurfTemp, iterNo, vaporOne, intermediateSaveOption, method); 

if compTwoVapMoF_s > 0
    ScCompTwo = calSchmidt(muVapFilm, rhoVapFilm, binDiffusivityTwo, ambient.Temperature, dropSurfTemp, iterNo, vaporTwo, intermediateSaveOption, method);
else
    ScCompTwo = 0 ; 
end

%% Sherwood number

Sh0CompOne = calSherwood(Re, Gr, ScCompOne, ambient.Temperature, dropSurfTemp, iterNo, vaporOne, intermediateSaveOption, sherwoodCorrelation);

if compTwoVapMoF_s > 0 
    Sh0CompTwo = calSherwood(Re, Gr, ScCompTwo, ambient.Temperature, dropSurfTemp, iterNo, vaporTwo, intermediateSaveOption, sherwoodCorrelation);
else
    Sh0CompTwo = 0 ; 
end

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if mikhil == 1 % Suitable for single and bi-component liquids !
    
    %% Calculate mass transfer spalding numbers
    
    vapMF_s = compOneVapMF_s + compTwoVapMF_s ; % surface vapor mass fraction
    
    BmCompOne = calSpaldingMT(compOneVapMF_s, vapMF_s, ambient.ComponentOneMF, ambient.Temperature, vaporOne, iterNo, intermediateSaveOption, method); % Spalding mass transfer number for water vapor

    if compTwoVapMoF_s > 0 
        BmCompTwo = calSpaldingMT(compTwoVapMF_s,vapMF_s, ambient.ComponentTwoMF, ambient.Temperature, vaporTwo, iterNo, intermediateSaveOption, method); % Spalding mass transfer number for urea vapor
    else
        BmCompTwo = 0 ; 
    end

    %% Modified Sherwood number

    ShStarCompOne = modifySherwood(Sh0CompOne, BmCompOne, Re, ambient.Temperature, dropSurfTemp, vaporOne, iterNo, intermediateSaveOption, method);

    if compTwoVapMoF_s > 0 
        ShStarCompTwo = modifySherwood(Sh0CompTwo, BmCompTwo, Re, ambient.Temperature, dropSurfTemp, vaporTwo, iterNo, intermediateSaveOption, method);
    else
        ShStarCompTwo = 0 ;
    end
    
    %% Equivalent Bm for the mixture (Sazhin et al. 2010)
    
    vapMF_f = ambient.ComponentOneMF + ambient.ComponentTwoMF ; % far field vapor mass fraction
    
    BmMix = calSpaldingMT(vapMF_s, vapMF_s, vapMF_f, ambient.Temperature, vaporMix, iterNo, intermediateSaveOption, method);
    
    %% Equivalent Diffusivity for the mixture (Sazhin et al. 2010)
    
    mixDiffusivity = calMixtureDiffusivity(ambient.Pressure, ambient.Temperature, dropSurfTemp, 0.333, compOneVapMoF_s, compTwoVapMoF_s, vaporOne, vaporTwo, ambient.Gas, method, diffusivityModel) ;
    
    %% Equivalent Schmidt number for the mixture
    
    ScMix = calSchmidt(muVapFilm, rhoVapFilm, mixDiffusivity, ambient.Temperature, dropSurfTemp, iterNo, vaporMix, intermediateSaveOption, method); 
    
    %% Equivalent Sh for the mixture
    
    Sh0Mix = calSherwood(Re, Gr, ScMix, ambient.Temperature, dropSurfTemp, iterNo, vaporMix, intermediateSaveOption, sherwoodCorrelation);
    
    ShStarMix = modifySherwood(Sh0Mix, BmMix, Re, ambient.Temperature, dropSurfTemp, vaporMix, iterNo, intermediateSaveOption, method);
    
    %% Equivalent Le for the mixture 
    
    LeMix = calLewis(condVapFilm, condAmbGas, cpFilm, cpAmbGas, rhoCompOneVap, rhoCompTwoVap, rhoAir, rhoVapFilm, mixDiffusivity, ambient.Temperature, dropSurfTemp, iterNo, vaporMix, intermediateSaveOption, method); 

    %% Calculate Heat Transfer Spalding Numbers
    
    if strcmp(componentTwoDepletionCalculationMethod,'usingVaporPressure') == 1 || (strcmp(componentTwoDepletionCalculationMethod,'usingArrheniusEquation') == 1 && dropSurfTemp < compOneBoilingPoint ) %le(compTwoLiqMaFr, 0.99999)) 
        
        [Bt, phi, NuStaar] = calSpaldingHT_iterative(cpVapor, cpFilm, cpAmbGas, ShStarMix, Nu0, Re, LeMix, BmMix, ambient.Temperature, dropSurfTemp, iterNo, intermediateSaveOption, vaporMix, method);
        
    else
        
        [Bt, phi, NuStaar] = calSpaldingHT_arrhenius(cpCompOneVap, cpFilm, mixEnthalpy, ShStarMix, Nu0, Re, LeMix, ambient.Temperature, dropSurfTemp, compOneBoilingPoint, iterNo, intermediateSaveOption, vaporOne, method);
        
        
    end
    
end