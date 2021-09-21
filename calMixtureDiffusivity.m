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

%% This function calculates an equivalent diffusivity of the binary vapour mixture into the ambient.
%  It may be suitably modified to use other expressions. 
%  Equivalent diffusivity is obtained using the method proposed by the work on bi-component droplet evaporation by Sazhin et al. (2010)

function [mixtureDiffusivity] = calMixtureDiffusivity(ambientPressure, ambientTemperature, surfaceTemperature, filmRulePara, compOneMoleFrac_s, compTwoMoleFrac_s, vapor1, vapor2, ambientGas, calMethod, diffModel)

%% select the method
mik = strcmp(calMethod, 'mikhil');

%% select the correlation
mikhil = strcmp(diffModel, 'mikhil');        % in effect there are only three options for the calculation of mixture diffusivity.  
fuller = strcmp(diffModel, 'fuller');
wilkeLee = strcmp(diffModel, 'wilkeLee');
chapmanEnskog = strcmp(diffModel, 'chapmanEnskog');

%% select the fluids
nitrogen = strcmp(ambientGas,'N2');
air = strcmp(ambientGas,'air');

waterVapor = strcmp(vapor1,'waterVapor');
rankWater = 1 ; 
if waterVapor ==0
    waterVapor = strcmp(vapor2,'waterVapor');
    rankWater = 2 ;
end

ureaVapor = strcmp(vapor1,'ureaVapor');
rankUrea = 1 ;
if ureaVapor ==0
    ureaVapor = strcmp(vapor2,'ureaVapor');
    rankUrea = 2 ; 
end

nHeptaneVapor = strcmp(vapor1,'nHeptaneVapor');
rankHeptane = 1 ;
if nHeptaneVapor ==0
    nHeptaneVapor = strcmp(vapor2,'nHeptaneVapor');
    rankHeptane = 2 ;
end

nDecaneVapor = strcmp(vapor1,'nDecaneVapor');
rankDecane = 1; 
if nDecaneVapor ==0
    nDecaneVapor = strcmp(vapor2,'nDecaneVapor');
    rankDecane = 2 ; 
end

%% calculate the binary diffusivity

if mik == 1 
    
    if wilkeLee == 1 
        
        filmTemperature = surfaceTemperature*(1-filmRulePara) + ambientTemperature*filmRulePara;
    
        a = 1.06036;
        b = 0.15610;
        c = 0.19300;
        d = 0.47635;
        e = 1.03587;
        f = 1.52996;
        g = 1.76474;
        h = 3.89411; 

        if nitrogen == 1   % Stephen Turns
    
            gas.MW = 28.014;
 
            gas.sigma = 3.798;
    
            gas.eByKB = 71.4;
    
        elseif air == 1    % Stephen Turns
    
            gas.MW = 28.97;
    
            gas.sigma = 3.711;
    
            gas.eByKB = 78.6;
    
        end

        if waterVapor == 1     % Stephen Turns
        
            if rankWater ==  1
                i = 1 ;
            else
                i = 2 ;
            end
    
            Vap(i).MW = 18.02;            
    
            Vap(i).sigma = 2.641;
    
            Vap(i).eByKB = 809.1;
   
    
        end

        if nHeptaneVapor == 1   % Estimation of LennardJones (6,12) Pair Potential Parameters from Gas Solubility Data - Wilhelm (1971)
        
            if rankHeptane ==  1
                i = 1 ;
            else
                i = 2 ;
            end
    
            Vap(i).MW = 100.204;            
    
            Vap(i).sigma = 6.25;
    
            Vap(i).eByKB = 573;
    
        end
    
        if nDecaneVapor == 1    % Estimation of LennardJones (6,12) Pair Potential Parameters from Gas Solubility Data - Wilhelm (1971)
        
            if rankDecane ==  1
                i = 1 ;
            else
                i = 2 ;
            end
    
            Vap(i).MW = 142.285;            
    
            Vap(i).sigma = 7.08;
    
            Vap(i).eByKB = 601;
    
        end

        if ureaVapor == 1     % Correlations used by Nagaraju
        
            if rankUrea ==  1
                i = 1 ;
            else
                i = 2 ;
            end
    
            Vap(i).MW = 60.06; 
    
            Vc = 218;
            Vb = 0.285*Vc^1.048;
            Vap(i).sigma = 1.18*Vb^(1/3);
    
            Tb = 465;
            Vap(i).eByKB = 1.15 * Tb;

        end
    
        weightOne = compOneMoleFrac_s / (compOneMoleFrac_s + compTwoMoleFrac_s) ;
        weightTwo = compTwoMoleFrac_s / (compOneMoleFrac_s + compTwoMoleFrac_s) ;
    
        vap.MW = weightOne * Vap(1).MW + weightTwo * Vap(2).MW ; 
    
        mwAB = 2/(1/vap.MW + 1/gas.MW);
    
        vap.sigma = weightOne * Vap(1).sigma + weightTwo * Vap(2).sigma ;  

        sigmaGV = (vap.sigma + gas.sigma)/2;
    
        eV = weightOne * Vap(1).eByKB + weightTwo * Vap(2).eByKB ;
    
        eG = gas.eByKB ; 
    
        eGV = (eV * eG)^0.5;
    
        TStar = filmTemperature / eGV; 
    
        omega = a/(TStar^b) + c/exp(d*TStar) + e/exp(f*TStar) + g/exp(h*TStar);
    
        mixtureDiffusivity = (( 3.03 - (0.98/mwAB^0.5))*10^(-7) * filmTemperature^1.5) / (ambientPressure*10^(-5) * mwAB^0.5 * sigmaGV^2 * omega ) ;  % Wilke and Lee Correlation
        
    end

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    if fuller == 1 || mikhil == 1        % using correlation by Fuller et al. Obtained from Perry's Handbook 2-370 and confirmed using the original paper. 
    
        filmTemperature = surfaceTemperature*(1-filmRulePara)+ ambientTemperature*filmRulePara;   

        if nitrogen == 1   
    
            gas.MW = 28.01;
 
            gas.nu = 17.9;
    
        end

        if air == 1  
    
            gas.MW = 28.97 ;
    
            gas.nu = 20.1 ;
    
        end

        if waterVapor == 1   
        
            if rankWater == 1 
                i = 1 ;
            else 
                i = 2 ; 
            end
    
            Vap(i).MW = 18.02;            
    
            Vap(i).nu = 12.7;
    
        end

        if nHeptaneVapor == 1  
        
            if rankHeptane == 1 
                i = 1 ;
            else 
                i = 2 ; 
            end
    
            Vap(i).MW = 100.204;            
    
            Vap(i).nu = 7 * 16.5 + 16 * 1.98 ;   % C7 H16
        
        end

        if nDecaneVapor == 1   
        
            if rankDecane == 1 
                i = 1 ;
            else 
                i = 2 ; 
            end
    
            Vap(i).MW = 142.285 ;            
    
            Vap(i).nu = 10 * 16.5 + 22 * 1.98 ;  % C10 H22 
    
        end

        if ureaVapor == 1 
        
            if rankUrea == 1 
                i = 1 ;
            else 
                i = 2 ; 
            end
    
            Vap(i).MW = 60.06; 
    
            Vap(i).nu = 16.5 + 4 * 1.98 + 2 * 5.69 + 5.481 ;  % C H4 N2 O 

        end
    
        weightOne = compOneMoleFrac_s / (compOneMoleFrac_s + compTwoMoleFrac_s) ;
        weightTwo = compTwoMoleFrac_s / (compOneMoleFrac_s + compTwoMoleFrac_s) ;
    
        vap.MW = weightOne * Vap(1).MW + weightTwo * Vap(2).MW ; 
    
        vap.nu = weightOne * Vap(1).nu + weightTwo * Vap(2).nu ; 
    
        mixtureDiffusivity = (0.01013 * filmTemperature^(1.75) * (1/vap.MW + 1/gas.MW)^(0.5))/(ambientPressure * ((vap.nu)^(1/3)+(gas.nu)^(1/3))^2) ; 

    end
    
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    if chapmanEnskog == 1 
        
        filmTemperature = surfaceTemperature*(1-filmRulePara)+ ambientTemperature*filmRulePara;

        a = 1.06036;
        b = 0.15610;
        c = 0.19300;
        d = 0.47635;
        e = 1.03587;
        f = 1.52996;
        g = 1.76474;
        h = 3.89411;

        if nitrogen == 1   % Stephen Turns
    
            gas.MW = 28.014;
 
            gas.sigma = 3.798;
    
            gas.eByKB = 71.4;

        elseif air == 1       % Stephen Turns
    
            gas.MW = 28.97;
    
            gas.sigma = 3.711;
    
            gas.eByKB = 78.6; 
    
        end

        if waterVapor == 1     % Stephen Turns
        
            if rankWater == 1 
                i = 1 ;
            else 
                i = 2 ; 
            end
    
            Vap(i).MW = 18.02;            
    
            Vap(i).sigma = 2.641;
    
            Vap(i).eByKB = 809.1;
    
        end

        if nHeptaneVapor == 1   % Estimation of LennardJones (6,12) Pair Potential Parameters from Gas Solubility Data - Wilhelm (1971)
        
            if rankHeptane == 1 
                i = 1 ;
            else 
                i = 2 ; 
            end
    
            Vap(i).MW = 100.204;            
    
            Vap(i).sigma = 6.25;
    
            Vap(i).eByKB = 573;
    
        end

        if nDecaneVapor == 1    % Estimation of LennardJones (6,12) Pair Potential Parameters from Gas Solubility Data - Wilhelm (1971)
        
            if rankDecane == 1 
                i = 1 ;
            else 
                i = 2 ; 
            end
    
            Vap(i).MW = 142.285;            
    
            Vap(i).sigma = 7.08;
    
            Vap(i).eByKB = 601;
    
        end

        if ureaVapor == 1     % Correlations used by Nagaraju
        
            if rankUrea == 1 
                i = 1 ;
            else 
                i = 2 ; 
            end
    
            Vap(i).MW = 60.06; 
    
            Vc = 218;
            Vb = 0.285*Vc^1.048;
            Vap(i).sigma = 1.18*Vb^(1/3);
    
            Tb = 465;
            Vap(i).eByKB = 1.15 * Tb;

        end
        
        weightOne = compOneMoleFrac_s / (compOneMoleFrac_s + compTwoMoleFrac_s) ;
        weightTwo = compTwoMoleFrac_s / (compOneMoleFrac_s + compTwoMoleFrac_s) ;
    
        vap.MW = weightOne * Vap(1).MW + weightTwo * Vap(2).MW ; 
        
        mwAB = 2/(1/vap.MW + 1/gas.MW);
        
        vap.sigma = weightOne * Vap(1).sigma + weightTwo * Vap(2).sigma ;

        sigmaGV = (vap.sigma + gas.sigma)/2;
        
        eV = weightOne * Vap(1).eByKB + weightTwo * Vap(2).eByKB ;
    
        eG = gas.eByKB ; 
    
        eGV = (eV * eG)^0.5;
    
        TStar = filmTemperature / eGV;
    
        omega = a/(TStar^b) + c/exp(d*TStar) + e/exp(f*TStar) + g/exp(h*TStar);
    
        mixtureDiffusivity = (0.0266 * filmTemperature^(3/2))/(ambientPressure * mwAB^(0.5) * sigmaGV^2 * omega) ;
        
    end
    
end



