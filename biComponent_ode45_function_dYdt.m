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

%% This function defines the derivatives of the variables that are being solved (for) by ODE45. 

function dy = biComponent_ode45_function_dYdt(t,y)

%% Read user inputs from file

inputFile ;           % read the user inputs from file

boilingPoints ;       % read the boiling points

constVelocity = strcmp(isVelocityConstant, 'yes');     % checks if the relative velocity between the droplet and the gas is constant throughout its lifetime. 

iterationNumber ;     % read the iteration number from file

usingVaporPressure = strcmp(componentTwoDepletionCalculationMethod,'usingVaporPressure');       % checks if the depletion of component two is to be modeled as an evaporation process
usingArrheniusEquation = strcmp(componentTwoDepletionCalculationMethod,'usingArrheniusEquation');       % checks if the depletion of component two is to be modeled using an Arrhenius equation model

%% Defining the variables

dSurfTemp = y(1) ;   % K

dia = y(2) ;  % m

volDrop = (1/6) * pi * dia^3 ; % m^3

areaDrop = pi * dia^2 ;         % surface area in m^2

compTwoLiqMF = y(3) ;

if compTwoLiqMF < 1e-5  % round off the liquid mass fraction of componente two. 
    
    compTwoLiqMF = 0 ; 
    
elseif compTwoLiqMF > 0.99999
    
    compTwoLiqMF = 1 ; 
    
end

relVelocity = y(4) ; % m/s

%% Evaluate and retrieve the properties that are needed to define the derivatives

[rhoLOne, dRhoOnedT, rhoLTwo, dRhoTwodT, rhoL, rhoFilm, binDiffOne, binDiffTwo, ShStarOne, ShStarTwo, BmOne, BmTwo, Bt, cpVap, condVap, mixEnth, cpL, Re, rhoInf, NuZero, NuStarr] = evaluateProperties(dSurfTemp, dia, compTwoLiqMF, relVelocity);
 
%% Defining other variables

massDrop = volDrop * rhoL ;

compTwoMass = massDrop * (compTwoLiqMF) ; 

%% Define the rate of change of mass, and the rate of change of diameter of the droplet

% rate of mass loss by evaporation of component one
if compTwoLiqMF < 1 
    
    if usingVaporPressure == 1  % This method models the evaporation of water as well as the depletion of urea as evaporation processess using their respective saturation pressures. 
        
        mDotOne = pi * dia * (rhoFilm * binDiffOne) * ShStarOne * log(1 + BmOne) ; % kg/s
        
    elseif usingArrheniusEquation == 1  % This method models the urea depletion process using an Arrhenius equation. Water evaporation is modeled using the saturation pressure of water. 
        
        if dSurfTemp < compOneBoilingPoint
            
            mDotOne = pi * dia * (rhoFilm * binDiffOne) * ShStarOne * log(1 + BmOne) ; % kg/s
            
        else
            
            mDotOne = pi * dia * (condVap/cpVap) * NuStarr(1,1) * log(1 + Bt) ;  % kg/s % Equation given by Birkhold et al. 10.1016/j.apcatb.2005.12.035
            % properties evaluated assuming that the droplet temperature remains constant once it reaches the boiling point of water (component one)
            
        end
        
    end
    
else
    mDotOne = 0 ;
end

% rate of mass loss by evaporation of component two
if compTwoLiqMF > 0 
    
    if usingVaporPressure == 1  % Models the evaporation of urea (or the second component) using its saturation pressure. 
        
        mDotTwo = pi * dia * (rhoFilm * binDiffTwo) * ShStarTwo * log(1 + BmTwo) ; % kg/s
        
    elseif usingArrheniusEquation == 1 % Depletion of urea modeled using an Arrhenius equation. (Arrhenius Equation can only be used for urea )
        
        if compTwoLiqMF < 1  % Urea depletion starts only after water evaporation is complete.  
            
            mDotTwo = 0 ; 
            
        else
        
            if strcmp(componentTwo.Name,'urea') == 1 
            
                [preExpFactor, actEnergy] = getArrheniusParameters(arrheniusParameterReference);  % arrheniusParameterReference is a user input. 
        
                mDotTwo = pi * dia * preExpFactor * exp(-actEnergy/(8.314 * dSurfTemp)) ;  % kg/s
            
            else 
            
                mDotTwo = 0 ;
            
            end
        
        end
        
    end
    
else
    mDotTwo = 0;
end

% effective rate of mass loss from the droplet 
dmdt = mDotOne + mDotTwo ;     

%% Rate of change of droplet temperature 

if usingVaporPressure == 1 

    if Bt > 0
    
        qDrop = dmdt * ( cpVap * ((ambient.Temperature - dSurfTemp)/Bt) - mixEnth ) ;  % An averaged Cp of the vapors is used here.
    
    elseif dSurfTemp < ambient.Temperature
    
        qDrop = ( NuZero * condVap / dia ) * (areaDrop) * (ambient.Temperature - dSurfTemp);  % To account for droplet heating when there is no evaporation. 
    
    else
    
        qDrop = 0 ; 
    
    end

end

% Birkhold et al. (10.1016/j.apcatb.2005.12.035) uses the following equation for evaluating qDrop to the urea particle
% Urea dissociation occurs only after the complete evaporation of water. 

if usingArrheniusEquation == 1  % (Arrhenius Equation can only be used for urea )

    if compTwoLiqMF < 1

        qDrop = dmdt * ( cpVap * ((ambient.Temperature - dSurfTemp)/Bt) - mixEnth ) ;  % An averaged Cp of the vapors is used here. (Same as using vapor pressure)

    else
    
        if dSurfTemp < ambient.Temperature
        
            if dmdt == 0 
            
                qDrop = ( NuZero * condVap / dia ) * (areaDrop) * (ambient.Temperature - dSurfTemp);  % To account for droplet heating when there is no evaporation.
            
            else
        
                qDrop = ( NuStarr(2,1) * condVap / dia ) * (areaDrop) * (ambient.Temperature - dSurfTemp) - dmdt * mixEnth ;  % kg/s % Equation given by Birkhold et al. 10.1016/j.apcatb.2005.12.035
            
            end
        
        else 
        
            qDrop = 0 ;
        
        end

    end
    
end


if strcmp(limitDropletTemperature, 'yes') == 1 % checks if the maximum droplet temperature is limited to the boiling point of either of the liquids. 
    
    if compTwoLiqMF < 1 && ge(dSurfTemp, min(compOneBoilingPoint, compTwoBoilingPoint))  % this condition may be met by UWS droplets evaporating at temperatures above BP of water. 
    
        dTdt = 0 ;
    
    else
    
        dTdt = qDrop / (volDrop * rhoL * cpL) ;
    
    end
    
else % if the droplet temperature is not limited to either of the boiling points. 
    
    dTdt = qDrop / (volDrop * rhoL * cpL) ;
    
end

%% Rate of change of the liquid mass fraction of component two

if compTwoLiqMF < 1 
    
    dYdt = (massDrop * (-mDotTwo) - compTwoMass * (-dmdt))/(massDrop^2) ;     % Using quotient rule for derivatives. dmdt and mDotTwo are multiplied with -1 because they are losses in mass
    
else 
    
    dYdt = 0 ; 
    
end

%% rate of change of liquid densities
%  The mixture density is a function of both temperature and mass fraction.
%  The rate of change of density is hence evaluated using the rate of change of density with temperature, and the rate of change of temperature and mass fraciton with respect to time.

dRhoOneBydTime = dRhoOnedT * dTdt ; 
dRhoTwoBydTime = dRhoTwodT * dTdt ; 

% calculate dRhodt
termOne = ( (1-compTwoLiqMF)/rhoLOne^2 ) * dRhoOneBydTime ;
if compTwoLiqMF > 0 
    termTwo = ( compTwoLiqMF/rhoLTwo^2 ) * dRhoTwoBydTime ;
    termThree = ( 1/rhoLOne - 1/rhoLTwo ) * (dYdt) ;
else
    termTwo = 0 ; 
    termThree = 0 ; 
end

dRhodt = (rhoL)^2 * (termOne + termTwo + termThree ) ;

%% rate of change of volume by time 

dVdt = ( rhoL * (-dmdt) - massDrop *  dRhodt )/ (rhoL)^2 ;              % using the quotient rule for derivatives. negative sign is added because it is a loss in mass         

dDdt = (6/pi)^(1/3) * (1/3) * (volDrop)^(-2/3) * dVdt ;                 % rate of change of diameter in m/s


%% Rate of change of droplet relative velocity 
%  Rate of change of velocity is evaluated when the velocity is not declared as a constant.  

if constVelocity ~= 1 && relVelocity > 1e-5

    Cd = calDragCoefficient(Re, ambient.Temperature, iterNo, intermediateSaveOption, method); 
    
    dUdt = ((3 * Cd) / (4 * dia)) * (rhoInf/rhoL) * relVelocity ^ 2 ; 
    
else 
    
    dUdt = 0 ;
    
end


%% Write iteration number to file

iterNo = iterNo + 1 ;

writeIterationNumber((iterNo));
    
dy = [dTdt; dDdt; dYdt; -dUdt]; % the derivatives of the variables that are being solved (for) by ODE45. 

