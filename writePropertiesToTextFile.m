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

%% This function evaluates and writes the values of various liquid and vapor thermophysical properties, as well as non-dimensional numbers.
%  The properties and non-dimensional numbers, at each time step, are evaluated and can be written to text files.  
%  This function is called after the solution is complete, and the data generated is one of the outputs of the code.

function writePropertiesToTextFile (time, temperature, diameter, compTwoLiqMassFraction, relativeVelocity)

inputFile ; 

boilingPoints ; 

numberOfEntries = size(temperature);

for iteration = 1 : numberOfEntries (1,1)
    
    %% Initialize the variables 
    
    dSurfTemp = temperature(iteration,1) ; 
    compTwoLiqMF = compTwoLiqMassFraction(iteration,1); 
    relVelocity = relativeVelocity(iteration,1);
    dia = diameter(iteration,1)/1000;
    
    %% Evaluate the properties 
    
    [compOneLiqMoFr, compTwoLiqMoFr, rhoLiqOne, dRhoOneBydTemp, rhoLiqTwo, dRhoTwoBydTemp, rhoLiq, cpLiqOne, cpLiqTwo, cpLiq, compOneSatPress, compTwoSatPress, totalSatPress, compOneVapMoF_s, compTwoVapMoF_s, vapMoF_s, ambGasMoF_s, compOneVapMF_s, compTwoVapMF_s, ambGasMF_s, compOneVapMoF_inf, compTwoVapMoF_inf, ambGasMoF_inf, compOneVapMF_film, compTwoVapMF_film, ambGasMF_film,  compOneVapMoF_film, compTwoVapMoF_film, ambGasMoF_film, filmMixMW, cpN2, cpCompOneVap, cpCompTwoVap, cpVapMix, condOneVap, condTwoVap, condN2, condO2, condAr, condAir, condAmbGas, condVapMix, diffOneInAir, diffTwoInAir, diffTwoInOne, muN2FilmTemp, muO2FilmTemp, muArFilmTemp, muAirFilmTemp, muAmbGas, muVapOne, muVapTwo, muN2AmbientTemp, muVapMix, rhoCompOneVap, rhoCompTwoVap, rhoInf, rhoAir, rhoVapMix, compOneEnth, compTwoEnth, mixEnth, cpVapOneTwo, LeCompOne, LeCompTwo, Re, Pr, ScCompOne, ScCompTwo, Sh0_CompOne, Sh0_CompTwo, Bm_CompOne, Bm_CompTwo, ShStar_CompOne, ShStar_CompTwo, Nu0, BT, PHI, NuSTAR] = evaluateProperties_postSolvingTheODE(dSurfTemp, dia, compTwoLiqMF, relVelocity);
       
    %% Write properties to file
   
    % liquid properties
    filename = 'output_liquidProperties_withHeader.txt';
    
    if iteration == 1 
        fileID = fopen(filename,'w');
        fprintf(fileID, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'Time(sec.)', 'moleFractionOne', 'moleFractionTwo','rhoOne(kg/m3)','dRhoOnedT(kg/m3 K)','rhoTwo(kg/m3)','dRhoTwodT(kg/m3 K)','rhoMix(kg/m3)', 'cpOne(J/mol.K)','cpTwo(J/mol.K)', 'cpMix(J/mol.K)', 'pSatOne(Pa)', 'pSatTwo(Pa)', 'pSatTotal(Pa)');
    else
        fileID = fopen(filename,'a');
    end
    
    fprintf(fileID, '%f\t%f\t%f\t%4.2f\t%f\t%4.2f\t%f\t%4.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\n', [time(iteration,1); compOneLiqMoFr; compTwoLiqMoFr; rhoLiqOne; dRhoOneBydTemp; rhoLiqTwo; dRhoTwoBydTemp; rhoLiq; cpLiqOne; cpLiqTwo; cpLiq; compOneSatPress; compTwoSatPress; totalSatPress]);

    fclose(fileID);
 
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % vapor properties
    filename = 'output_vaporProperties_withHeader.txt';
    
    if iteration == 1 
        fileID = fopen(filename,'w');
        fprintf(fileID, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'Time(sec.)', 'MoleFracOne(s)','MoleFracTwo(s)', 'MoleFracVap(s)', 'MoleFracGas(s)', 'MassFracOne(s)', 'MassFracTwo(s)', 'MassFracGas(s)', 'MoleFracOne(inf)','MoleFracTwo(inf)', 'MoleFracGas(inf)', 'MassFracOne(film)','MassFracTwo(film)', 'MassFracGas(film)', 'MoleFracOne(film)','MoleFracTwo(film)', 'MoleFracGas(film)', 'mixMolWt', 'cpGas(J/kg.K)', 'cpOne(J/kg.K)', 'cpTwo(J/kg.K)', 'cpVapOneTwo(J/kg.K)','cpFilm(J/kg.K)', 'condOne(W/m.K)', 'condTwo(W/m.K)', 'condN2(W/m.K)', 'condO2(W/m.K)', 'condAr(W/m.K)', 'condAir(W/m.K)', 'condAmbGas(W/m.K)', 'condMix(W/m.K)', 'diffOneInAir(m2/s)', 'diffTwoInAir(m2/s)', 'diffTwoInOne(m2/s)', 'muN2FilmTemp(Pa.s)', 'muO2FilmTemp(Pa.s)', 'muArFilmTemp(Pa.s)', 'muAirFilmTemp(Pa.s)', 'muAmbGas(Pa.s)', 'muVapOne(Pa.s)', 'muVapTwo(Pa.s)', 'muN2AmbientTemp(Pa.s)', 'muVapMix(Pa.s)', 'rhoOne(kg/m3)', 'rhoTwo(kg/m3)', 'rhoInf(kg/m3)', 'rhoAir(kg/m3)', 'rhoFilm(kg/m3)', 'enthalpyOne(J/mol)', 'enthalpyTwo(J/mol)', 'enthalpyMix(J/kg)');
    else
        fileID = fopen(filename,'a');
        fprintf(fileID, '%f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%3.9f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%9.3f\t%9.3f\t%9.3f\n', [time(iteration,1); compOneVapMoF_s; compTwoVapMoF_s; vapMoF_s; ambGasMoF_s; compOneVapMF_s; compTwoVapMF_s; ambGasMF_s; compOneVapMoF_inf; compTwoVapMoF_inf; ambGasMoF_inf; compOneVapMF_film; compTwoVapMF_film; ambGasMF_film;  compOneVapMoF_film; compTwoVapMoF_film; ambGasMoF_film; filmMixMW; cpN2; cpCompOneVap; cpCompTwoVap; cpVapOneTwo; cpVapMix; condOneVap; condTwoVap; condN2;  condO2; condAr; condAir; condAmbGas; condVapMix; diffOneInAir; diffTwoInAir; diffTwoInOne; muN2FilmTemp; muO2FilmTemp; muArFilmTemp; muAirFilmTemp; muAmbGas; muVapOne; muVapTwo; muN2AmbientTemp; muVapMix; rhoCompOneVap; rhoCompTwoVap; rhoInf; rhoAir; rhoVapMix; compOneEnth; compTwoEnth; mixEnth]);
    end
    
    fclose(fileID);
    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
    % non dimensional numbers

    Nu_BP = NuSTAR(1,1) ;
    Nu_Ts = NuSTAR(2,1) ;
    
    filename = 'output_nonDimensionalNumbers_withHeader.txt';

    if iteration == 1 
        fileID = fopen(filename,'w');
        fprintf(fileID, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'Time(sec.)', 'LeOne','LeTwo', 'rel.Vel(m/s)', 'Re', 'Pr', 'ScOne', 'ScTwo', 'Sh0_one','Sh0_Two', 'BmOne', 'BmTwo', 'Sh*One', 'Sh*Two', 'Nu0', 'Bt', 'phi', 'Nu*(BP)', 'Nu*(Ts)');
    else
        fileID = fopen(filename,'a');
    end
    
    fprintf(fileID, '%f\t%6.6f\t%6.6f\t%6.6f\t%9.3f\t%8.4f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%9.3f\t%9.3f\t%6.6f\t%6.6f\n', [time(iteration,1); LeCompOne; LeCompTwo; relVelocity; Re; Pr; ScCompOne; ScCompTwo; Sh0_CompOne; Sh0_CompTwo; Bm_CompOne; Bm_CompTwo; ShStar_CompOne; ShStar_CompTwo; Nu0; BT; PHI; Nu_BP; Nu_Ts]);

    fclose(fileID);

end


