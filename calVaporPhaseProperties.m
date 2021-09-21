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

%% This function evaluates the vapor phase thermophysical properties based on the ambient conditions.
%  The properties are returned to the main function and can be saved to files.   

function [compOneVapMoleFrac_s, compTwoVapMoleFrac_s, vapMoleFrac_s, ambGasMoleFrac_s, compOneVapMassFrac_s, compTwoVapMassFrac_s, ambGasMassFrac_s, compOneVapMoleFrac_inf, compTwoVapMoleFrac_inf, ambGasMoleFrac_inf, compOneVapMassFrac_film, compTwoVapMassFrac_film, ambientGasMassFrac_film, compOneVapMoleFrac_film, compTwoVapMoleFrac_film, ambGasMoleFrac_film, filmMixMolWt, cpAmbient, cpCompOneVapor, cpCompTwoVapor, cpVapFilm, condOneVapor, condTwoVapor, condNitrogen, condOxygen, condArgon, cond_Air, condAmbient, condVaporFilm, binDiffusivityVapOne, binDiffusivityVapTwo, binDiffVapTwoVapOne, muNitrogen, muOxygen, muArgon, muAir, muAmbient, muVaporOne, muVaporTwo, muN2AmbientTemp, muVaporFilm, rhoCompOneVapor, rhoCompTwoVapor, rhoAirInf, rhoAirAtFilmTemp, rhoVaporFilm, compOneEnthalpy, compTwoEnthalpy, evapMixtureEnthalpy] = calVaporPhaseProperties(ambPressure, ambTemperature, dropSurfaceTemp, saturationPressureCompOne, saturationPressureCompTwo, liqCompTwoMassFrac, compOneAmbientMassFrac, compTwoAmbientMassFrac, compOneMolWt, compTwoMolWt, ambientGas, ambGasMolWt, vapOne, vapTwo, gasLaw, I, writeData, calculationMethod, binDiffModel, isEnthDissociationIncluded, ureaVaporizationEnthalpyCorrln, compTwoDepletionCalculationMethod)

FRP_vapor = 0.333;
% Liquid phase properties are evaluated using the one third rule.

% Liquid phase properties are evaluated at the droplet (surface) temperature
FRP_liquid = 0;

% Ambient gas properties are evaluated at the ambient temperature. 
FRP_gas = 1; 

boilingPoints ; 

saveDataToFile = strcmp(writeData,'yes');   % checks if intermediate data needs to be saved to files after each iteration. 

savePath = './Properties/vaporPhase/';

%% Vapor Mole Fractions at the droplet surface

compOneVapMoleFrac_s = calVapMoleFrac(ambPressure, saturationPressureCompOne, calculationMethod) ; 

if liqCompTwoMassFrac > 0 
    compTwoVapMoleFrac_s = calVapMoleFrac(ambPressure, saturationPressureCompTwo, calculationMethod) ;
else 
    compTwoVapMoleFrac_s = 0 ;
end

vapMoleFrac_s = compOneVapMoleFrac_s + compTwoVapMoleFrac_s;  % total mole fraction of the vapors

ambGasMoleFrac_s = 1 - (vapMoleFrac_s) ;   % mole fraction of the ambient gas

% write data to files
if saveDataToFile == 1
    
    filename = [savePath,'compOneMolFrac','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Mol Frac');
    end    
    fprintf(fileID, '%d\t %f\n', I, compOneVapMoleFrac_s);
    fclose(fileID);

    filename = [savePath,'compTwoMolFrac','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Mol Frac');
    end    
    fprintf(fileID, '%d\t %f\n', I, compTwoVapMoleFrac_s);
    fclose(fileID);

    filename = [savePath,'vapMolFrac','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Mol Frac');
    end    
    fprintf(fileID, '%d\t %f\n', I, vapMoleFrac_s);
    fclose(fileID);

    filename = [savePath,'ambGasMolFrac','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Mol Frac');
    end    
    fprintf(fileID, '%d\t %f\n', I, ambGasMoleFrac_s);
    fclose(fileID);
    
end

%% Vapor Mass Fractions at the droplet surface

[compOneVapMassFrac_s, compTwoVapMassFrac_s] = calSurfaceVapMassFrac(compOneVapMoleFrac_s, compTwoVapMoleFrac_s, compOneMolWt, compTwoMolWt, ambGasMolWt, calculationMethod) ;

ambGasMassFrac_s = 1 - (compOneVapMassFrac_s + compTwoVapMassFrac_s);

% write data to files
if saveDataToFile == 1
    
    filename = [savePath,'compOneVapMassFrac','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Mass Frac');
    end    
    fprintf(fileID, '%d\t %f\n', I, compOneVapMassFrac_s);
    fclose(fileID);

    filename = [savePath,'compTwoVapMassFrac','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Mass Frac');
    end    
    fprintf(fileID, '%d\t %f\n', I, compTwoVapMassFrac_s);
    fclose(fileID);

    filename = [savePath,'ambGasMassFrac','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Mass Frac');
    end    
    fprintf(fileID, '%d\t %f\n', I, ambGasMassFrac_s);
    fclose(fileID);
    
end

%% Vapor mole fractions far from the droplet 

[compOneVapMoleFrac_inf, compTwoVapMoleFrac_inf, ambGasMoleFrac_inf] = calVapMoleFrac_ambient(compOneAmbientMassFrac, compTwoAmbientMassFrac, compOneMolWt, compTwoMolWt, ambGasMolWt) ;

%% Average mass fractions in the film

[compOneVapMassFrac_film, compTwoVapMassFrac_film, ambientGasMassFrac_film] = calFilmMassFractions(compOneVapMassFrac_s, compTwoVapMassFrac_s, compOneAmbientMassFrac, compTwoAmbientMassFrac, calculationMethod);  


%% Average mole fractions in the film 

[compOneVapMoleFrac_film, compTwoVapMoleFrac_film, ambGasMoleFrac_film] = calVapMoleFrac_film(compOneVapMassFrac_film, compTwoVapMassFrac_film, compOneMolWt, compTwoMolWt, ambGasMolWt);


%% Molecular weight of mixture in the film

filmMixMolWt = calFilmMixtureMolecularWeight(compOneVapMoleFrac_film, compTwoVapMoleFrac_film, ambGasMoleFrac_film, compOneMolWt, compTwoMolWt, ambGasMolWt, calculationMethod); 

% write data to files
if saveDataToFile == 1
    
    filename = [savePath,'filmMixMolWt','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'MW (g/mol)');
    end    
    fprintf(fileID, '%d\t %f\n', I, filmMixMolWt);
    fclose(fileID);
    
end

%% Vapor phase specific heats

cpCompOneVapor = calVaporSpecificHeat(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, vapOne, calculationMethod); % J/kg K

if liqCompTwoMassFrac > 0 
    cpCompTwoVapor = calVaporSpecificHeat(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, vapTwo, calculationMethod); % J/kg K
else 
    cpCompTwoVapor = 0 ;
end

cpAmbient = calVaporSpecificHeat(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, ambientGas, calculationMethod); % J/kg K  

cpVapFilm = calVaporFilmSpecificHeat(cpCompOneVapor, cpCompTwoVapor, cpAmbient, compOneVapMassFrac_s, compOneAmbientMassFrac, compTwoVapMassFrac_s, compTwoAmbientMassFrac, FRP_vapor, calculationMethod); % J / kg K

% write data to files
if saveDataToFile == 1
    
    filename = [savePath,'cpVaporOne','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Cp (J/kg K)');
    end    
    fprintf(fileID, '%d\t %f\n', I, cpCompOneVapor);
    fclose(fileID);

    filename = [savePath,'cpVaporTwo','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Cp (J/kg K)');
    end    
    fprintf(fileID, '%d\t %f\n', I, cpCompTwoVapor);
    fclose(fileID);

    filename = [savePath,'cpNitrogen','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');   
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Cp (J/kg K)');
    end    
    fprintf(fileID, '%d\t %f\n', I, cpAmbient);
    fclose(fileID);

    filename = [savePath,'cpVaporFilm','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Cp (J/kg K)');
    end    
    fprintf(fileID, '%d\t %f\n', I, cpVapFilm);
    fclose(fileID);
    
end

%% Vapor phase viscosities 

muNitrogen = calViscosity(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, 'N2', calculationMethod);
muOxygen = calViscosity(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, 'O2', calculationMethod);
muArgon = calViscosity(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, 'Ar', calculationMethod);

moleFracN2 = 0.78 ; % mole fraction of nitrogen in air
moleFracO2 = 0.21 ; % mole fraction of oxygen in air
moleFracAr = 0.01 ; % mole fraction of argon in air
molarMassO2 = 32 ;  % molecular weight of oxygen
molarMassAr = 39.95; % molecular weight of argon

muAir = calFilmViscosity(muArgon, muOxygen, muNitrogen, moleFracAr, moleFracO2, moleFracN2, molarMassAr, molarMassO2, 'N2', calculationMethod); % viscosity of air is calculated assuming it to be a mixture of nitrogen, oxygen and argon. 

muVaporOne = calViscosity(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, vapOne, calculationMethod); 

if liqCompTwoMassFrac > 0 
    muVaporTwo = calViscosity(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, vapTwo, calculationMethod);
else 
    muVaporTwo = 0 ;
end

muN2AmbientTemp = calViscosity(ambPressure, ambTemperature, dropSurfaceTemp, FRP_gas,'N2', calculationMethod);

if strcmp(ambientGas,'N2')== 1       % The ambient gas may be N2 or air. 
    muAmbient = muNitrogen ;   
else
    muAmbient = muAir ;   
end

% viscosity of the vapor film. 
muVaporFilm = calFilmViscosity(muVaporOne, muVaporTwo, muAmbient, compOneVapMoleFrac_film, compTwoVapMoleFrac_film, ambGasMoleFrac_film, compOneMolWt, compTwoMolWt, ambientGas, calculationMethod);

% write data to files
if saveDataToFile == 1
    
    filename = [savePath,'muVaporOne','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'mu (Pa s)');
    end    
    fprintf(fileID, '%d\t %f\n', I, muVaporOne);
    fclose(fileID);

    filename = [savePath,'muVaporTwo','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'mu (Pa s)');
    end    
    fprintf(fileID, '%d\t %f\n', I, muVaporTwo);
    fclose(fileID);

    filename = [savePath,'muNitrogen','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'mu (Pa s)');
    end    
    fprintf(fileID, '%d\t %f\n', I, muNitrogen);
    fclose(fileID);

    filename = [savePath,'muN2Reynolds','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'mu (Pa s)');
    end    
    fprintf(fileID, '%d\t %f\n', I, muN2AmbientTemp);
    fclose(fileID);

    filename = [savePath,'muVaporFilm','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'mu (Pa s)');
    end    
    fprintf(fileID, '%d\t %f\n', I, muVaporFilm);
    fclose(fileID);
    
end




%% Vapor phase conductivities

condOneVapor = calVaporConductivity(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, vapOne, calculationMethod);
condNitrogen = calVaporConductivity(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, 'N2', calculationMethod);
condOxygen = calVaporConductivity(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, 'O2', calculationMethod);
condArgon = calVaporConductivity(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, 'O2', calculationMethod);

% conductivity of air as a mixture of nitrogen, oxygen and argon. 
cond_Air = calFilmConductivity(condArgon, condOxygen, condNitrogen, muArgon, muOxygen, muNitrogen, moleFracAr, moleFracO2, moleFracN2, molarMassAr, molarMassO2, 'N2', calculationMethod); 

if liqCompTwoMassFrac > 0 
    condTwoVapor = calVaporConductivity(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, vapTwo, calculationMethod);
else 
    condTwoVapor = 0 ;
end

if strcmp(ambientGas,'N2')== 1       % The ambient gas may be N2 or air. 
    condAmbient = condNitrogen ;   
else
    condAmbient = cond_Air ;   
end

% conductivity of the vapor film. 
condVaporFilm = calFilmConductivity(condOneVapor, condTwoVapor, condAmbient, muVaporOne, muVaporTwo, muNitrogen, compOneVapMoleFrac_film, compTwoVapMoleFrac_film, ambGasMoleFrac_film, compOneMolWt, compTwoMolWt, ambientGas, calculationMethod); 
   

% write data to files
if saveDataToFile == 1

    filename = [savePath,'condVaporOne','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'K (W/m K)');
    end    
    fprintf(fileID, '%d\t %f\n', I, condOneVapor);
    fclose(fileID);

    filename = [savePath,'condVaporTwo','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'K (W/m K)');
    end    
    fprintf(fileID, '%d\t %f\n', I, condTwoVapor);
    fclose(fileID);

    filename = [savePath,'condNitrogen','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'K (W/m K)');
    end    
    fprintf(fileID, '%d\t %f\n', I, condNitrogen);
    fclose(fileID);

    filename = [savePath,'condVaporFilm','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'K (W/m K)');
    end    
    fprintf(fileID, '%d\t %f\n', I, condVaporFilm);
    fclose(fileID);
    
end



%% Binary diffusivity

binDiffusivityVapOne = calBinaryDiffusivity(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, vapOne, 'N2', binDiffModel);

if liqCompTwoMassFrac > 0 
    binDiffusivityVapTwo = calBinaryDiffusivity(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, vapTwo, 'N2', binDiffModel);
    binDiffVapTwoVapOne = calBinaryDiffusivity(ambPressure, ambTemperature, dropSurfaceTemp, FRP_vapor, vapTwo, vapOne, binDiffModel);
else 
    binDiffusivityVapTwo = 0 ;
    binDiffVapTwoVapOne = 0;
end

% write data to files
if saveDataToFile == 1
    
    filename = [savePath,'diffVaporOneNitrogen','_',num2str(ambTemperature),'K.txt'];  
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Dab (m^2/s)');
    end    
    fprintf(fileID, '%d\t %f\n', I, binDiffusivityVapOne);
    fclose(fileID);

    filename = [savePath,'diffVaporTwoNitrogen','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Dab (m^2/s)');
    end    
    fprintf(fileID, '%d\t %f\n', I, binDiffusivityVapTwo);
    fclose(fileID);

    filename = [savePath,'diffVaporTwoVaporOne','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Dab (m^2/s)');
    end    
    fprintf(fileID, '%d\t %f\n', I, binDiffVapTwoVapOne);
    fclose(fileID);
    
end


% Corrected binary diffusivities (So as to account for the presence of the third component)

if compOneVapMassFrac_s > 0 && compTwoVapMassFrac_s > 0 % binary diffusivities are corrected so as to account for the presence of the third component. S. R. Turns (3rd Edition, page no. 228)
    
    L = zeros(3,3) ; 
    
    Chi = zeros(3,1) ; % film mole fraction matrix 
    Chi(1,1) = compOneVapMoleFrac_film ; 
    Chi(2,1) = compTwoVapMoleFrac_film ;
    Chi(3,1) = ambGasMoleFrac_film ;
    
    MW = zeros(3,1) ; % molecular weight matrix 
    MW(1,1) = compOneMolWt ; 
    MW(2,1) = compTwoMolWt ;
    MW(3,1) = ambGasMolWt ; 
    
    DAB = zeros(3,3) ; % diffusivity matrix
    DAB(2,1) = binDiffVapTwoVapOne; % diffusivity of the vapor of component two in the vapor of component one 
    DAB(1,3) = binDiffusivityVapOne; % diffusivity of the vapor of component one in the ambient gas
    DAB(2,3) = binDiffusivityVapTwo; % diffusivity of the vapor of component two in the ambient gas
    DAB(1,2) = DAB(2,1) ;
    DAB(3,1) = DAB(1,3) ;
    DAB(3,2) = DAB(2,3) ;
    
    for lc = 1:3  % find inverse of the correction matrix 
        for cm = 1:3
            
            if lc ~= cm
                
                lm = 6 - (lc + cm) ; 
                
                L(lc,cm) = Chi(cm,1) * ( MW(cm,1) * Chi(cm,1) + MW(lc,1) * Chi(lc,1) )/ ( MW(lc,1) * DAB (lc,cm) ) + Chi(lm,1) * MW(cm,1) * Chi(cm,1) / ( MW(lc,1) * DAB (lc,lm) ) ; 
                
            else
                
                L(lc,cm) = 0 ; 
                
            end
            
        end
    end
     
    isColumnZero = find(all(L==0)) ;
    
    if isColumnZero == 0 % ascertain that matrix "L" has no columns with all elements as zero  
        
        F = inv(L) ;  % find the correction matrix 
    
        for lc = 1:3 
            for cm = 1:3
            
                DAB(lc,cm) = Chi(lc,1) * filmMixMolWt / MW(cm,1) * (F(lc,cm) - F(lc,lc)); 
            
            end
        end
    
        binDiffusivityVapOne = DAB(1,3) ; % diffusivity of the vapor of component one in the ambient gas in the presence of component two
        binDiffusivityVapTwo = DAB(2,3) ; % diffusivity of the vapor of component two in the ambient gas in the presence of component one
        binDiffVapTwoVapOne  = DAB(2,1) ; % diffusivity of the vapor of component two in the vapor of component two in the presence of the ambient gas
    
    end
    
end



%% Vapor phase densities

rhoAirAtFilmTemp = calRhoGas(ambPressure,ambTemperature,dropSurfaceTemp,FRP_vapor, ambientGas, gasLaw, calculationMethod); % using ideal / real gas law
% Viscosities of other vapors may be obtained by changing the last input of
% the above function. Available options are air, N2, waterVapor, and
% ureaVapor.

rhoCompOneVapor = calRhoGas(ambPressure,ambTemperature,dropSurfaceTemp,FRP_vapor, vapOne, gasLaw, calculationMethod);

if liqCompTwoMassFrac > 0 
    rhoCompTwoVapor = calRhoGas(ambPressure,ambTemperature,dropSurfaceTemp,FRP_vapor, vapTwo, gasLaw, calculationMethod);
else 
    rhoCompTwoVapor = 0 ;
end

rhoAirInf = calRhoGas(ambPressure,ambTemperature,dropSurfaceTemp,FRP_gas, ambientGas, gasLaw, calculationMethod);
% viscosity and density of air for the evaluation of reynolds number is
% evaluated at the ambient gas temperature.

% density of the film. 
rhoVaporFilm = calRhoFilm(ambPressure, ambTemperature, dropSurfaceTemp, compOneMolWt, compTwoMolWt, compOneVapMoleFrac_s, compTwoVapMoleFrac_s, ambGasMoleFrac_s, rhoCompOneVapor, rhoCompTwoVapor, rhoAirAtFilmTemp, compOneVapMassFrac_film, compTwoVapMassFrac_film, ambientGasMassFrac_film, FRP_vapor, ambientGas, gasLaw, calculationMethod);

% write data to files
if saveDataToFile == 1
    
    filename = [savePath,'rhoVaporOne','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'rho (kg / m^3)');
    end    
    fprintf(fileID, '%d\t %f\n', I, rhoCompOneVapor);
    fclose(fileID);

    filename = [savePath,'rhoVaporTwo','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'rho (kg / m^3)');
    end    
    fprintf(fileID, '%d\t %f\n', I, rhoCompTwoVapor);
    fclose(fileID);

    filename = [savePath,'rhoAir','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'rho (kg / m^3)');
    end    
    fprintf(fileID, '%d\t %f\n', I, rhoAirAtFilmTemp);
    fclose(fileID);

    filename = [savePath,'rhoAirReynolds','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'rho (kg / m^3)');
    end    
    fprintf(fileID, '%d\t %f\n', I, rhoAirInf);
    fclose(fileID);

    filename = [savePath,'rhoVaporFilm','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'rho (kg / m^3)');
    end    
    fprintf(fileID, '%d\t %f\n', I, rhoVaporFilm);
    fclose(fileID);
    
end

%% Enthalpy of vaporization 

compOneEnthalpy = calVaporizationEnthalpy(ambTemperature, dropSurfaceTemp, FRP_liquid, vapOne, compOneBoilingPoint, calculationMethod, isEnthDissociationIncluded, ureaVaporizationEnthalpyCorrln);  % J/mol

if liqCompTwoMassFrac > 0 
    compTwoEnthalpy = calVaporizationEnthalpy(ambTemperature, dropSurfaceTemp, FRP_liquid, vapTwo, compTwoBoilingPoint, calculationMethod, isEnthDissociationIncluded, ureaVaporizationEnthalpyCorrln); % J/mol 
else 
    compTwoEnthalpy = 0 ;
end

if strcmp(compTwoDepletionCalculationMethod,'usingVaporPressure') == 1
    
    evapMixtureEnthalpy = calEnthalpyOfEvaporatingVapor(compOneEnthalpy, compTwoEnthalpy, compOneVapMoleFrac_s, compTwoVapMoleFrac_s, compOneMolWt, compTwoMolWt, calculationMethod); % J/kg
    
elseif strcmp(compTwoDepletionCalculationMethod,'usingArrheniusEquation') == 1  % In this case urea evaporation starts only after the completion of evaporation of water
    
    if liqCompTwoMassFrac < 0.99999  
    
        evapMixtureEnthalpy = compOneEnthalpy / compOneMolWt * 1000 ; 
        
    else
        
        evapMixtureEnthalpy = compTwoEnthalpy / compTwoMolWt * 1000 ;
        
    end
    
end   

% write data to files
if saveDataToFile == 1
    
    filename = [savePath,'enthalpyCompOne','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Hvap (J / mol)');
    end    
    fprintf(fileID, '%d\t %f\n', I, compOneEnthalpy);
    fclose(fileID);

    filename = [savePath,'enthalpyCompTwo','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Hvap (J / mol)');
    end    
    fprintf(fileID, '%d\t %f\n', I, compTwoEnthalpy);
    fclose(fileID);

    filename = [savePath,'enthalpyEvapMixture','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Hvap (J / kg)');
    end    
    fprintf(fileID, '%d\t %f\n', I, evapMixtureEnthalpy);
    fclose(fileID);
    
end
