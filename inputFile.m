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

%% --------------------User-specified inputs----------------------- %%
droplet.InitialDia = 0.001;    % in m
droplet.InitialSurfaceTemp = 300;    % in K

componentOne.Name = 'water';    % options are 'water', 'nHeptane', 'nDecane'
componentTwo.Name ='urea';    % options are 'urea', 'water', 'nHeptane', 'nDecane' (component with lower volatility has to be component two)
componentOne.MW = 18.02 ;    % in g/mol water 18.02, urea 60.06, nHeptane 100.204, nDecane 142.285.
componentTwo.MW = 60.06 ;    % in g/mol

componentTwo.initMassFrac = 0.325 ; 
partialPressureLaw = 'raoults';    % options are 'raoults' and 'NRTL'. NRTL can be used only for UWS. 

ambient.Temperature = 673;    % in K  
ambient.Pressure = 101325;    % in Pa 
ambient.Gas = 'N2' ;    % N2 and air are the options.  
ambient.GasMW = 28.01 ;    % in g/mol
ambient.ComponentOneMF = 0 ;    % value between 0 and 1 
ambient.ComponentTwoMF = 0 ;    % value between 0 and 1 
ambient.airRelativeVelocity = 4.3 ;    % initial velocity in m/s

isVelocityConstant = 'yes' ;    % if the relative velocity between air and droplet is constant, then set this value to 'yes' else set it to 'no'. 

equationOfState = 'ideal' ;    % options are 'ideal' and 'real'. 

liquidPropertyCharacter = 'temperatureDependent';    % options are 'temperatureDependent' and 'temperatureIndependent' (properties are either evaluated at a constant mean temperature or at actual droplet temperature)  

diffusivityModel = 'fuller';    % options are 'chapmanEnskog', 'fuller', 'wilkeLee' and 'mikhil'.  

sherwoodCorrelation = 'ranzMarshall';    % options are 'ranzMarshall', 'renksizbulut', 'frossling', and 'mikhil'.

nusseltCorrelation  = 'ranzMarshall';    % options are 'ranzMarshall', 'renksizbulut', 'frossling', and 'mikhil'.

isNaturalConvectionIncluded = 'no' ;    % options are 'yes' or 'no'. (applicable when the relative velocity between air and the droplet becomes zero)

componentTwoDepletionCalculationMethod = 'usingVaporPressure' ; % options are 'usingVaporPressure' and 'usingArrheniusEquation'

limitDropletTemperature = 'no' ;  % setting this to yes limits the maximum droplet temperature to the boiling point of the component with higher volatility

startTime = 0 ;       % start of computation time. always set this to zero 
endTime   = 400 ;    % end of computation time in seconds  

intermediateSaveOption = 'no';    % set to 'yes' to write intermediate property values to files (saving will make the code very slow and was included only for debugging)

method ='mikhil';    %  the method is adapted from the method proposed by Abramzon and Sirignano 

%% -------------------------UWS Specific Inputs--------------------------------------------------------------------------------------------------------------------------------------------

ureaVaporPressureCorrelation = 'ebrahimian' ; % options are 'bernhard', birkhold, 'ebrahimian', 'ferrq' and 'lundstrom'  

ureaVaporizationEnthalpyCorrelation = 'lundstrom' ;    % options are 'birkhold' (temperature independent) and 'lundstrom' (temperature dependent)
isDissociationEnergyIncluded = 'no' ;    % options are 'yes' or 'no'. (applicable only for UWS) 

uwsDensityCalculationMethod = 'moltenUrea' ;     % options are 'constant', 'moltenUrea'.   
uwsSpHeatCalculationMethod = 'moltenUrea' ;    % option currently available is 'moltenUrea'.

arrheniusParameterReference = 'birkhold' ;     % option currently available is 'birkhold'

%-------------------------------------------------------------------------%

