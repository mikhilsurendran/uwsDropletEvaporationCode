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

%% This function modifies the sherwood number.  
%  Can be modified to use other expressions for modifying Sh
%  Inputs and outputs are in SI units

function [ShStaar] = modifySherwood(ShO, Bm, Reynolds, ambientTemperature, surfacefTemperature, vapor, I, writeData, calMethod)

saveDataToFile = strcmp(writeData,'yes');
savePath = './Properties/nonDimen/';

%% select the correlation
mikhil = strcmp(calMethod, 'mikhil');

%% select the fluid
waterVapor = strcmp(vapor,'waterVapor');
ureaVapor = strcmp(vapor,'ureaVapor');
nHeptaneVapor = strcmp(vapor,'nHeptaneVapor');
nDecaneVapor = strcmp(vapor,'nDecaneVapor');
vaporMixture = strcmp(vapor,'vaporMix');

%% calculate the modified Sherwood number

if mikhil == 1 
    
    if ureaVapor == 1
    
        if ge(surfacefTemperature, 406.15)
        
            if Bm > 0 
            
                fBM = (((1+Bm)^0.7) * log(1+Bm)) / Bm ;
            
                ShStaar = 2 + (ShO - 2) / fBM ;
            
            else 
            
                ShStaar = ShO ; 
            
            end
        
        else
        
            ShStaar = 0 ; 
        
        end
    
        % write data to file
        if saveDataToFile == 1
        
            filename = [savePath,'modShUrea','_',num2str(ambientTemperature),'K.txt'];
            fileID = fopen(filename,'a');
            if I == 1
                fprintf(fileID, '%s\t %s\n', 'i', 'modified Sh');
            end    
            fprintf(fileID, '%d\t %f\n', I, ShStaar);
            fclose(fileID);
        
        end

    else
    
        if Bm > 0 
        
            fBM = (((1+Bm)^0.7) * log(1+Bm)) / Bm ;
          
            ShStaar = 2 + (ShO - 2) / fBM ;
          
        else
        
            ShStaar = ShO;
        
        end
    
        % write data to file
        if saveDataToFile == 1
        
            if waterVapor == 1
                filename = [savePath,'modShWater','_',num2str(ambientTemperature),'K.txt'];
            
            elseif nHeptaneVapor == 1
                filename = [savePath,'modShHeptane','_',num2str(ambientTemperature),'K.txt'];
            
            elseif nDecaneVapor == 1
                filename = [savePath,'modShDecane','_',num2str(ambientTemperature),'K.txt'];
                
            elseif vaporMixture == 1
                filename = [savePath,'modShMixture','_',num2str(ambientTemperature),'K.txt'];
            
            end
               
            fileID = fopen(filename,'a');
            if I == 1
                fprintf(fileID, '%s\t %s\n', 'i', 'modified Sh');
            end    
            fprintf(fileID, '%d\t %f\n', I, ShStaar);
            fclose(fileID);
        
        end
    
    end
        
end