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

%% This function writes output data to text files (after the solution is complete)

function writeOutputDataToTextFile (time, temperature, diameter, dropletMass, componentTwoMass, componentTwoMassFraction,  deltaM, deltaMOne, deltaMTwo, relativeVelocity)

inputFile ; 

time = time';
temperature = temperature';
diameter = diameter';
dropletMass = dropletMass' ;
componentTwoMass = componentTwoMass' ;
componentTwoMassFraction = componentTwoMassFraction';
relativeVelocity = relativeVelocity';

d2Byd02 = diameter.^2 / (diameter(1)^2) ;

componentOneMass = dropletMass - componentTwoMass ; 

deltaM = deltaM' ; 
deltaMOne = deltaMOne' ; 
deltaMTwo = deltaMTwo' ; 


% output data with header
filename = 'outputFile_dropletData_withHeader.txt';

fileID = fopen(filename,'w');

fprintf(fileID, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'Time (s)', 'Temperature (K)','Diameter(mm)', 'd^2/do^2', 'Droplet Mass (ug)', 'Mass of component one (ug)', 'Mass of component two (ug)', 'Component two Mass Fraction', 'dMassDt (ug/s)', 'dMOneDt (ug/s)', 'dMTwoDt (ug/s)', 'Droplet Relative Velocity(m/s)');
fprintf(fileID, '%4.8f\t%4.3f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.11f\t%1.4f\t%1.11f\t%1.11f\t%1.11f\t%2.6f\n', [time; temperature; diameter; d2Byd02; dropletMass; componentOneMass; componentTwoMass; componentTwoMassFraction; deltaM; deltaMOne; deltaMTwo; relativeVelocity]);

fclose(fileID);

%% user-specified inputs are also written to a file. 

initDia = droplet.InitialDia ;
initTemp = droplet.InitialSurfaceTemp ; 
nameOne = componentOne.Name' ;
nameTwo = componentTwo.Name ;
compOneMW = componentOne.MW ;
compTwoMW = componentTwo.MW ;
ambTemp = ambient.Temperature ;
ambPress = ambient.Pressure ;
ambGas = ambient.Gas ;

filename = 'outputFile_inputConditions.txt';

fileID = fopen(filename,'w');

fprintf(fileID, '%s\t%1.8f\n', 'Initial Diameter (m) - ', initDia);
fprintf(fileID, '%s\t%4.3f\n', 'Initial Drop Temperature (K) - ', initTemp);
fprintf(fileID, '%s\t%s\n', 'Component One - ', nameOne');
fprintf(fileID, '%s\t%s\n', 'Component Two - ', nameTwo');
fprintf(fileID, '%s\t%3.3f\n', 'Mol Weight One - ', compOneMW);
fprintf(fileID, '%s\t%3.3f\n', 'Mol Weight Two - ', compTwoMW);
fprintf(fileID, '%s\t%1.4f\n', 'Component TWo Mass Fraction - ', componentTwo.initMassFrac);
fprintf(fileID, '%s\t%s\n', 'Component Two Depletion Method - ', componentTwoDepletionCalculationMethod');
fprintf(fileID, '%s\t%s\n', 'Urea Vapor Pressure Correlation - ', ureaVaporPressureCorrelation');
fprintf(fileID, '%s\t%s\n', 'Urea Vaporization Enthalpy Correlation - ', ureaVaporizationEnthalpyCorrelation');
fprintf(fileID, '%s\t%s\n', 'Is Dissociation Energy Included? - ', isDissociationEnergyIncluded');
fprintf(fileID, '%s\t%s\n', 'Partial Pressure Law - ', partialPressureLaw');
fprintf(fileID, '%s\t%4.3f\n', 'Ambient Temperature (K) - ', ambTemp);
fprintf(fileID, '%s\t%9.3f\n', 'Ambient Pressure (Pa) - ', ambPress);
fprintf(fileID, '%s\t%s\n', 'Ambient Gas - ', ambGas');
fprintf(fileID, '%s\t%1.3f\n', 'Initial Rel. Vel. (m/s) - ', ambient.airRelativeVelocity);
fprintf(fileID, '%s\t%s\n', 'Is Rel. Vel. Constant? - ', isVelocityConstant');
fprintf(fileID, '%s\t%s\n', 'Diffusivity Model - ', diffusivityModel');
fprintf(fileID, '%s\t%s\n', 'Sh Corrln - ', sherwoodCorrelation');
fprintf(fileID, '%s\t%s\n', 'Nu Corrln - ', nusseltCorrelation');
fprintf(fileID, '%s\t%s\n', 'Is Natural Convection Included? - ', isNaturalConvectionIncluded');
fprintf(fileID, '%s\t%s\n', 'uws Density - ', uwsDensityCalculationMethod');
fprintf(fileID, '%s\t%s\n', 'uws Sp. Heat - ', uwsSpHeatCalculationMethod');
fprintf(fileID, '%s\t%s\n', 'Liquid Properties - ', liquidPropertyCharacter');

fclose(fileID);