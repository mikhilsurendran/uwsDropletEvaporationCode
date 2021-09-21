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

%% This function calculates the enthalpy of vaporization. It may be suitably modified to use other expressions. 

function [enthalpy] = calVaporizationEnthalpy(ambientTemperature, surfaceTemperature, filmRulePara, compound, boilPoint, calMethod, isEnthDissociationAdded, ureaVapEnthalpyCorrln)

%% select the correlation
mikhil = strcmp(calMethod, 'mikhil');

%% select the correlation for the enthalpy of vaporization for urea
birkhold = strcmp (ureaVapEnthalpyCorrln, 'birkhold');
lundstrom = strcmp(ureaVapEnthalpyCorrln, 'lundstrom');

%% select the fluid
water = strcmp(compound,'waterVapor');
urea = strcmp(compound,'ureaVapor');
nHeptane = strcmp(compound,'nHeptaneVapor');
nDecane = strcmp(compound,'nDecaneVapor');

%% check whether the enthalpy of dissociation is to be added to the enthalpy of vaporization 

addDissociation = strcmp (isEnthDissociationAdded, 'yes'); 

%% Calculate the latent heat of vaporation based on selected correlation

if mikhil == 1 

    filmTemperature = surfaceTemperature*(1-filmRulePara)+ ambientTemperature*filmRulePara;  % filmRulePara is set to zero. So the vaporization enthalpies of the liquids are evaluated at the liquid temperature (surfaceTemperature). 

    if water == 1        % Correlation From Perry's Handbook, not from original paper.
    
        Tc = getCriticalTemperature('waterVapor') ;
        
        if filmTemperature > boilPoint 
            
            filmTemperature = boilPoint ; 
            
        end
    
        Tr = filmTemperature / Tc ;
    
        C1 = 5.2053 * 10^7; 
        C2 = 0.3199 ;
        C3 = -0.212 ;
        C4 = 0.25795 ;
    
        enthalpy = C1 * (1-Tr)^(C2 + C3*Tr + C4*Tr^2); % J/kmol
    
        enthalpy = enthalpy /1000 ; % J/mol
    

    elseif urea == 1
        
        if birkhold == 1   % used by Birkhold et al. as well as Nishad et al. 
                
            if addDissociation == 1
            
                enthalpy = ( 8.74 + 9.81 )* 10^4 ; % J/mol 
                
            else
                
                enthalpy = ( 8.74 )* 10^4 ; % J/mol 
                
            end
                
        elseif lundstrom == 1
                
            if addDissociation == 1   % as described in the paper by Lundstrom in DOI: 10.1177/0954407011406048
                    
                A = 24.856;
                B = 1.4437*10^(-1);
                C = 3.8088*10^(-5);
                D = -1.1007*10^(-7);
                E = 3.9161*10^(-11);
                
                meltingPoint = 406 ; 
             
                constantOne = 78210 ; % enthalpy of vaporization of urea at its melting point
                
                % Cp_ureaVapor = A + (B) * (filmTemperature) + (C) * (filmTemperature^2) + (D) * (filmTemperature^3) + (E) * (filmTemperature^4);    %J/mol K  % correlation given in Chemical Properties Handbook by Carl L. Yaws 
                
                % integrating the above equation between 0 K and filmTemperature, we have 

                % integralCpdT_ureaVapor = A * (filmTemperature) + (B/2) * (filmTemperature^2) + (C/3) * (filmTemperature^3) + (D/4) * (filmTemperature^4) + (E/5) * (filmTemperature^5);   %J/mol K
                
                % substituting the melting point of urea, we have 
                
                constantTwo = A * (meltingPoint) + (B/2) * (meltingPoint^2) + (C/3) * (meltingPoint^3) + (D/4) * (meltingPoint^4) + (E/5) * (meltingPoint^5) ;  % J/mol K                 
                
                a = 965.507;
                b = -5.0993;
                c = 1.0028*10^(-2);
                d = -6.3799*10^(-6); 
                
                % Cp_ureaLiquid = a + (b) * (filmTemperature) + (c) * (filmTemperature^2) + (d) * (filmTemperature^3) ;   %J/mol K  % correlation given in Chemical Properties Handbook by Carl L. Yaws 
                
                % integrating the above equation between 0 K and filmTemperature, we have 
                
                % integralCpdT_ureaLiquid = a * (filmTemperature) + (b/2) * (filmTemperature^2) + (c/3) * (filmTemperature^3) + (d/4) * (filmTemperature^4)   %J/mol K 
                
                % substituting the melting point of urea, we have 
                
                constantThree = a * (meltingPoint) + (b/2) * (meltingPoint^2) + (c/3) * (meltingPoint^3) + (d/4) * (meltingPoint^4) ;  %J/mol K 
                
                constantFour = constantOne - constantTwo + constantThree ;
                
                enthalpy = constantFour + ( (A-a) * (filmTemperature) + (B-b)/2 * (filmTemperature^2) + (C-c)/3 * (filmTemperature^3) + (D-d)/4 * (filmTemperature^4) + (E/5) * (filmTemperature^5) ) ; % J/mol 
            
                enthalpy = enthalpy + 9.81 * 10^4 ; % J/mol 
                
            else
                    
                A = 24.856;
                B = 1.4437*10^(-1);
                C = 3.8088*10^(-5);
                D = -1.1007*10^(-7);
                E = 3.9161*10^(-11);
                
                meltingPoint = 406 ; 
             
                constantOne = 78210 ; % enthalpy of vaporization of urea at its melting point
                
                % Cp_ureaVapor = A + (B) * (filmTemperature) + (C) * (filmTemperature^2) + (D) * (filmTemperature^3) + (E) * (filmTemperature^4);    %J/mol K  % correlation given in Chemical Properties Handbook by Carl L. Yaws 
                
                % integrating the above equation between 0 K and filmTemperature, we have 

                % integralCpdT_ureaVapor = A * (filmTemperature) + (B/2) * (filmTemperature^2) + (C/3) * (filmTemperature^3) + (D/4) * (filmTemperature^4) + (E/5) * (filmTemperature^5);   %J/mol K
                
                % substituting the melting point of urea, we have 
                
                constantTwo = A * (meltingPoint) + (B/2) * (meltingPoint^2) + (C/3) * (meltingPoint^3) + (D/4) * (meltingPoint^4) + (E/5) * (meltingPoint^5) ;  % J/mol K                 
                
                a = 965.507;
                b = -5.0993;
                c = 1.0028*10^(-2);
                d = -6.3799*10^(-6); 
                
                % Cp_ureaLiquid = a + (b) * (filmTemperature) + (c) * (filmTemperature^2) + (d) * (filmTemperature^3) ;   %J/mol K  % correlation given in Chemical Properties Handbook by Carl L. Yaws 
                
                % integrating the above equation between 0 K and filmTemperature, we have 
                
                % integralCpdT_ureaLiquid = a * (filmTemperature) + (b/2) * (filmTemperature^2) + (c/3) * (filmTemperature^3) + (d/4) * (filmTemperature^4)   %J/mol K 
                
                % substituting the melting point of urea, we have 
                
                constantThree = a * (meltingPoint) + (b/2) * (meltingPoint^2) + (c/3) * (meltingPoint^3) + (d/4) * (meltingPoint^4) ;  %J/mol K 
                
                constantFour = constantOne - constantTwo + constantThree ;
                
                enthalpy = constantFour + ( (A-a) * (filmTemperature) + (B-b)/2 * (filmTemperature^2) + (C-c)/3 * (filmTemperature^3) + (D-d)/4 * (filmTemperature^4) + (E/5) * (filmTemperature^5) ) ;
                
            end
                
        end
    
    
    elseif nHeptane == 1   % Correlation From Perry's Handbook, not from original paper.
    
        Tc = getCriticalTemperature('nHeptaneVapor') ;
        
        if filmTemperature > boilPoint 
            
            filmTemperature = boilPoint ; 
            
        end
    
        Tr = filmTemperature / Tc ;
    
        C1 = 5.0014 * 10^7; 
        C2 = 0.38795 ;
        C3 = 0 ;
        C4 = 0 ;
    
        enthalpy = C1 * (1-Tr)^(C2 + C3*Tr + C4*Tr^2); % J/kmol
    
        enthalpy = enthalpy /1000 ; % J/mol
    
    elseif nDecane == 1   % Correlation from original paper
        
        Tc = getCriticalTemperature('nDecaneVapor') ;
        
        if filmTemperature > boilPoint 
            
            filmTemperature = boilPoint ; 
            
        end
        
        molecularWt = 142.285 ; % kg/kmol
    
        C1 = 9.453 ; 
        C2 = 0.38 ;
    
        enthalpy = C1 * (Tc - filmTemperature)^C2 ; % kcal/kg
    
        enthalpy = ( enthalpy * 4.184 * molecularWt ); % J/mol
    
    end

end