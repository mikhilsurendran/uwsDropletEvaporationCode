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

%% This function evaluates the liquid phase thermophysical properties based on the ambient conditions.
%  The properties are returned to the main function and can be saved to files.   

function [rhoLiquidOne, dRhoOneBydTemp, rhoLiquidTwo, dRhoTwoBydTemp, rhoLiquid, cpLiquidOne, cpLiquidTwo, cpLiquid, satPressureCompOne, satPressureCompTwo, totalVaporPressure] = calLiquidPhaseProperties(ambPressure, ambTemperature, initialDropTemperature, dropSurfaceTemp, liqCompTwoMassFrac, liqCompTwoMoleFrac, liqCompOneMoleFrac, compOneMolWt, compTwoMolWt, I, writeData, compOneName, compTwoName, compOneBoilingPoint, compTwoBoilingPoint, vaporPressureLaw, calculationMethod, character, exptalOrCorrln, ureaPsatCorrln, compTwoDepletionCalculationMethod)

FRP_liquid = 0;
% Liquid phase properties are evaluated at the droplet (surface) temperature
% Vapor phase properties are evaluated using the one third rule. FRP_vapor = 0.333;
% Ambient gas properties are evaluated at the ambient temperature. FRP_gas = 1;

liqCompOneMassFrac = 1 - liqCompTwoMassFrac ; 

saveDataToFile = strcmp(writeData,'yes');  % checks if intermediate data needs to be saved to files after each iteration. 

savePath = './Properties/liquidPhase/';

%% Liquid phase densities

saltName = compTwoName ; % this variable is used to identify urea. 
saltMF = liqCompTwoMassFrac ; 

[rhoLiquidOne, dRhoOneBydTemp] = calRhoLiquid(ambPressure, ambTemperature, initialDropTemperature, compOneBoilingPoint, dropSurfaceTemp, FRP_liquid, compOneName, saltName, saltMF, calculationMethod, character, exptalOrCorrln);

if liqCompTwoMassFrac > 0  
    [rhoLiquidTwo, dRhoTwoBydTemp] = calRhoLiquid(ambPressure, ambTemperature, initialDropTemperature, compTwoBoilingPoint, dropSurfaceTemp, FRP_liquid, compTwoName, saltName, saltMF, calculationMethod, character, exptalOrCorrln);
else 
    rhoLiquidTwo = 0;
    dRhoTwoBydTemp = 0 ; 
end

% liquid mixture density
rhoLiquid = calMixDensity(rhoLiquidOne, rhoLiquidTwo, liqCompTwoMassFrac, calculationMethod);

if saveDataToFile == 1
    
    % write data to files
    
    filename = [savePath,'compOneDensity','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'rho (kg/m^3)');
    end    
    fprintf(fileID, '%d\t %f\n', I, rhoLiquidOne);
    fclose(fileID);

    filename = [savePath,'compTwoDensity','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'rho (kg/m^3)');
    end 
    fprintf(fileID, '%d\t %f\n', I, rhoLiquidTwo);
    fclose(fileID);
    
    filename = [savePath,'mixtureDensity','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'rho (kg/m^3)');
    end 
    fprintf(fileID, '%d\t %f\n', I, rhoLiquid);  
    fclose(fileID);
    
end

%% Liquid phase specific heats

cpLiquidOne = calLiquidSpecificHeat(ambPressure, ambTemperature, initialDropTemperature, compOneBoilingPoint, dropSurfaceTemp, FRP_liquid, compOneName, calculationMethod, character);

if liqCompTwoMassFrac > 0
    cpLiquidTwo = calLiquidSpecificHeat(ambPressure, ambTemperature, initialDropTemperature, compTwoBoilingPoint, dropSurfaceTemp, FRP_liquid, compTwoName, calculationMethod, character);
else
    cpLiquidTwo = 0;
end

% liquid mixture specific heat
cpLiquid = liqCompOneMassFrac*cpLiquidOne + liqCompTwoMassFrac*cpLiquidTwo;

if saveDataToFile == 1

    % write data to files
    filename = [savePath,'compOneSpHeat','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Cp (J/mol K)');
    end 
    fprintf(fileID, '%d\t %f\n', I, cpLiquidOne);
    fclose(fileID);

    filename = [savePath,'compTwoSpHeat','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Cp (J/mol K)');
    end 
    fprintf(fileID, '%d\t %f\n', I, cpLiquidTwo);
    fclose(fileID);

    filename = [savePath,'mixtureSpHeat','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Cp (J/mol K)');
    end 
    fprintf(fileID, '%d\t %f\n', I, cpLiquid);   
    fclose(fileID);
    
end

%% Vapor saturation pressures at the droplet surface


vapOne = strcat(compOneName,'Vapor');
vapTwo = strcat(compTwoName,'Vapor');

if liqCompTwoMassFrac < 1
    satPressureCompOne = calSatPressure(ambTemperature, dropSurfaceTemp, compOneBoilingPoint, FRP_liquid, vapOne, calculationMethod, ureaPsatCorrln) ; 
else 
    satPressureCompOne = 0;
end

if strcmp(compTwoDepletionCalculationMethod,'usingVaporPressure') == 1
    if liqCompTwoMassFrac > 0
        satPressureCompTwo = calSatPressure(ambTemperature, dropSurfaceTemp, compTwoBoilingPoint, FRP_liquid, vapTwo, calculationMethod, ureaPsatCorrln) ;
    else
        satPressureCompTwo = 0;
    end
elseif strcmp(compTwoDepletionCalculationMethod,'usingArrheniusEquation') == 1  % Arrhenius equation can be used only for urea
    satPressureCompTwo = 0;
end
    
% partial pressure
[satPressureCompOne, satPressureCompTwo, totalVaporPressure] = calPartialPressure(satPressureCompOne, satPressureCompTwo, ambTemperature, dropSurfaceTemp, liqCompTwoMassFrac, compOneBoilingPoint, FRP_liquid, compOneMolWt, compTwoMolWt, vaporPressureLaw, calculationMethod, compTwoName) ;


if saveDataToFile == 1

    % write data to files
    filename = [savePath,'compOnePsat','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Psat (Pa)');
    end 
    fprintf(fileID, '%d\t %f\n', I, satPressureCompOne);
    fclose(fileID);

    filename = [savePath,'compTwoPsat','_',num2str(ambTemperature),'K.txt'];
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Psat (Pa)');
    end 
    fprintf(fileID, '%d\t %f\n', I, satPressureCompTwo);
    fclose(fileID);
    
end

