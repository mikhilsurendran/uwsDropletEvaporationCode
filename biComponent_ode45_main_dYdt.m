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

%% The main file 

tic;  

clc ; 
clear all;  

%% Reading user inputs

inputFile; 

initiateCount = 1; 

writeIterationNumber(initiateCount);


%% Calculate the boiling points (used to calculate fluid properties by Abramzon and Sirignano)

[componentOne.BoilingPoint, componentTwo.BoilingPoint] = calBoilingPoint(componentOne.Name, componentTwo.Name, componentTwo.initMassFrac, ambient.Pressure);

lowerBP = min(componentOne.BoilingPoint, componentTwo.BoilingPoint) ;

writeBoilingPoint (componentOne.BoilingPoint, componentTwo.BoilingPoint);


%% Inputs for ODE45

initialConditions = [droplet.InitialSurfaceTemp droplet.InitialDia componentTwo.initMassFrac ambient.airRelativeVelocity]; 

timeSpan = [startTime endTime];


%% Solve the equations using ODE45
 
options = odeset('Events',@eventsFunction,'RelTol',1e-11,'AbsTol',[1e-11 1e-11 1e-11 1e-11]);

[t,y] = ode45(@biComponent_ode45_function_dYdt, timeSpan, initialConditions, options); % time (t) is the independent variable, and y (T, d, Yu, U) are the depedent variables. 

%% Print on screen

toc;

disp('*************************************************************************************************************************************************************');
disp('*** Solution Complete *** Solution Complete *** Solution Complete *** Solution Complete *** Solution Complete *** Solution Complete *** Solution Complete ***'); 
disp('*************************************************************************************************************************************************************');

% disp('Press any key to continue the program !')  
% pause;

%% Fetch the derivatives of the variables stored in y

dy = zeros(length(y),4) ; 

for count = 1:length(y)
    
    dy(count,:) = biComponent_ode45_function_dYdt(t(count,1), y(count,:));   % the derivates of T, d, Yu and U are stored in dy. 
    
end

%% conditioning the data and writing output data to a text file

droplet.Temperature = y(:,1);
droplet.Diameter = y(:,2) ;                         
droplet.ComponentTwoMF = y(:,3); 
droplet.D2 = droplet.Diameter .^2;
droplet.normalizedD2 = droplet.D2 ./ ((droplet.InitialDia)^2);
droplet.relativeVelocity = y(:,4);
droplet.dTempDt = dy(:,1);
droplet.dDiaDt = dy(:,2);
droplet.dY2Dt = dy(:,3);
droplet.dVelDt = dy(:,4); 

%calculate droplet mass and its rate of change
[droplet.Mass, droplet.ComponentOneMass, droplet.ComponentTwoMass, droplet.dMassDtime, droplet.dMassOneDtime, droplet.dMassTwoDtime] = calDropletMass(droplet.Diameter, droplet.Temperature, componentOne.BoilingPoint, componentTwo.BoilingPoint, droplet.ComponentTwoMF, componentOne.Name, componentTwo.Name, ambient.Pressure, ambient.Temperature, droplet.dTempDt, droplet.dDiaDt, droplet.dY2Dt, method, liquidPropertyCharacter, uwsDensityCalculationMethod);                

droplet.Diameter = droplet.Diameter * 1e3 ; % converting from m to mm 
droplet.D2 = droplet.D2 * 1e6 ;             % converting from m^2 to mm^2
normalizedTime = t ./ ((droplet.InitialDia * 1000 )^2) ; 

writeOutputDataToTextFile(t, droplet.Temperature, droplet.Diameter, droplet.Mass, droplet.ComponentTwoMass, droplet.ComponentTwoMF, droplet.dMassDtime, droplet.dMassOneDtime, droplet.dMassTwoDtime, droplet.relativeVelocity);
writePropertiesToTextFile(t, droplet.Temperature, droplet.Diameter, droplet.ComponentTwoMF, droplet.relativeVelocity); 

%% plotting data  

plotData(t, droplet.Diameter, droplet.InitialDia, droplet.ComponentTwoMF, droplet.Temperature, droplet.D2, droplet.normalizedD2, normalizedTime, droplet.Mass, droplet.ComponentOneMass, droplet.ComponentTwoMass, droplet.dMassDtime, droplet.dMassOneDtime, droplet.dMassTwoDtime, lowerBP, componentTwo.Name);

beep;
 