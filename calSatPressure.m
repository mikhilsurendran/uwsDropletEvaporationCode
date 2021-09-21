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

%% This function calculates the saturation pressure for the vapors using correlations. It may be suitably modified to use other expressions. 
%  By default droplet surface temperature is used to evaluate the
%  saturation pressure. Hence filmRulePara is set to zero. 

function [satPressure] = calSatPressure(ambientTemperature, surfaceTemperature, boilingPoint, filmRulePara, gas, calMethod, ureaCorrelation)

%% select the liquid
waterVapor = strcmp(gas,'waterVapor');
ureaVapor = strcmp(gas,'ureaVapor');
nHeptaneVapor = strcmp(gas,'nHeptaneVapor');
nDecaneVapor = strcmp(gas,'nDecaneVapor');

%% select the correlation
mikhil = strcmp (calMethod, 'mikhil');

%% select the correlation for the vapor pressure of urea
bernhard = strcmp(ureaCorrelation, 'bernhard');
birkhold = strcmp(ureaCorrelation, 'birkhold');
ebrahimian = strcmp(ureaCorrelation, 'ebrahimian');
ferrq = strcmp(ureaCorrelation, 'ferrq');
lundstrom = strcmp(ureaCorrelation, 'lundstrom');

%% set the temperature for property evaluation

if ureaVapor ~= 1      
    
    % this subroutine is applicable only for liquids. not applicable for salts like urea.
    if le(surfaceTemperature, boilingPoint)
        
        temperature = surfaceTemperature ; 
        
    else
        
        temperature = boilingPoint ;    % if the evaluated droplet temperature is higher than the boiling point of a component, then its properties are evaluated at its boiling point. 
        
    end

else
    
    temperature = surfaceTemperature ; 
    
end

%% Calculate the saturation vapor pressure based on the selected correlation

if mikhil == 1

    if waterVapor == 1         % Correlation From Perry's Handbook, not from original paper.
    
        C1 = 73.649;
        C2 = -7258.2;
        C3 =-7.3037;
        C4 = 4.1653e-06;
        C5 = 2;
        satPressure = exp(C1+C2/temperature+C3*log(temperature)+C4*temperature^C5);
  
    end

    if ureaVapor == 1
        
        if bernhard == 1 % used by Nishad et al.
            
            if ge(surfaceTemperature, 406.15)
            
                satPressure = exp(32.472-11755/temperature);
            
            else 
            
                satPressure = 0 ; 
            
            end
            
        elseif birkhold == 1 % used by Birkhold et al.
            
            if ge(surfaceTemperature, 406.15)
            
                satPressure = exp(12.06-3992/temperature);
            
            else 
            
                satPressure = 0 ; 
            
            end
            
        elseif ebrahimian == 1 % used by Ebrahimian et al. 
            
            if ge(surfaceTemperature, 406.15)
            
                satPressure = exp(29.9548-10876.1/temperature);
            
            else 
            
                satPressure = 0 ; 
            
            end
            
        elseif ferrq == 1 % as described in the paper by Lundstrom in DOI: 10.1177/0954407011406048
            
            if ge(surfaceTemperature, 406.15)
                
                A = 24.856;
                B = 1.4437*10^(-1);
                C = 3.8088*10^(-5);
                D = -1.1007*10^(-7);
                E = 3.9161*10^(-11);
                
                meltingPoint = 406 ; 
             
                constantOne = 78210 ; % enthalpy of vaporization of urea at its melting point
                
                % Cp_ureaVapor = A + (B) * (temperature) + (C) * (temperature^2) + (D) * (temperature^3) + (E) * (temperature^4);    %J/mol K % correlation given in Chemical Properties Handbook by Carl L. Yaws 
                
                % integrating the above equation between 0 K and T, we have 

                % integralCpdT_ureaVapor = A * (T) + (B/2) * (T^2) + (C/3) * (T^3) + (D/4) * (T^4) + (E/5) * (T^5);   %J/mol K
                
                % substituting the melting point of urea, we have 
                
                constantTwo = A * (meltingPoint) + (B/2) * (meltingPoint^2) + (C/3) * (meltingPoint^3) + (D/4) * (meltingPoint^4) + (E/5) * (meltingPoint^5) ;  % J/mol K                 
                
                a = 965.507;
                b = -5.0993;
                c = 1.0028*10^(-2);
                d = -6.3799*10^(-6); 
                
                % Cp_ureaLiquid = a + (b) * (temperature) + (c) * (temperature^2) + (d) * (temperature^3) ;   %J/mol K  % correlation given in Chemical Properties Handbook by Carl L. Yaws 
                
                % integrating the above equation between 0 K and T, we have 
                
                % integralCpdT_ureaLiquid = a * (T) + (b/2) * (T^2) + (c/3) * (T^3) + (d/4) * (T^4)   %J/mol K 
                
                % substituting the melting point of urea, we have 
                
                constantThree = a * (meltingPoint) + (b/2) * (meltingPoint^2) + (c/3) * (meltingPoint^3) + (d/4) * (meltingPoint^4) ;  %J/mol K 
                
                constantFour = constantOne - constantTwo + constantThree ; 
                
                % vaporizationEnthalpy = constantFour + (A * (T) + (B/2) * (T^2) + (C/3) * (T^3) + (D/4) * (T^4) + (E/5) * (T^5)) - ( a * (T) + (b/2) * (T^2) + (c/3) * (T^3) + (d/4) * (T^4)) ;
                
                % vaporizationEnthalpy = constantFour + ( (A-a) * (T) + (B-b)/2 * (T^2) + (C-c)/3 * (T^3) + (D-d)/4 * (T^4) + (E/5) * (T^5) ) ;
                
                R = 8.314 ; % J / mol K
                
                % integrating (vaporizationEnthalpy / RT^2) between 406 K and T we have, 
                
                % integralVapEnthalpy = (constantFour/R)*(1/406 - 1/T) + 1/R * {(A-a) * ln(T) + ((B-b)/2) * T + ((C-c)/6) * T^2 + ((D-d)/12) * T^3 + (E/20) * T^4 } - constantFive 
                
                % constantFive = 1/R * ( (A-a) * ln(meltingPoint) + ((B-b)/2) * meltingPoint + ((C-c)/6) * meltingPoint^2 + ((D-d)/12) * meltingPoint^3 + (E/20) * meltingPoint^4 )
                
                constantFive = 1/R * ( (A-a) * log(meltingPoint) + ((B-b)/2) * meltingPoint + ((C-c)/6) * meltingPoint^2 + ((D-d)/12) * meltingPoint^3 + (E/20) * meltingPoint^4 ) ;
                
                termOne = ( constantFour / R ) * (1/meltingPoint - 1/temperature) ;
                
                termTwo = 1/R * ( (A-a) * log(temperature) + ((B-b)/2) * temperature + ((C-c)/6) * temperature^2 + ((D-d)/12) * temperature^3 + (E/20) * temperature^4 ) ;
                
                satPressure = ( 10 ^ (10.3 - 4750 / meltingPoint) ) * 10^3 * exp( termOne + termTwo - constantFive) ; % in Pa
            
            else 
            
                satPressure = 0 ; 
            
            end
            
        elseif lundstrom == 1  % using the simplified correlation given by by Lundstrom in DOI: 10.1177/0954407011406048
            
            if ge(surfaceTemperature, 406.15)
                
                A = 11.663; 
                B = 4085.4; 
                C = 3.3953 * 10^-3; 
                
                satPressure = 10 ^ ( A - B / (C+temperature)) ;
                
            else
                
                satPressure = 0 ;
                
            end
            
        end
    
    end

    if nHeptaneVapor == 1  
    
        % Correlation from Abramzon's original paper (1989).
        C1 = 10.5892;
        C2 = -3934.731;
   
        satPressure = exp(C1+C2/temperature);
        satPressure = satPressure * 101325 ;
    
    end
    
    if nDecaneVapor == 1
    
        % Correlation from Abramzon's original paper (1989).
        C1 = 11.495;
        C2 = -5141.36;
   
        satPressure = exp(C1+C2/temperature);
        satPressure = satPressure * 101325 ; 

    end
    
end