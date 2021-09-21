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

%% This function iteratively calculates the Spalding heat transfer number.
%  Inputs and outputs are in SI units

function [Bt, phi, NuStar] = calSpaldingHT_iterative(CpV, CpFilm, CpG, ShStaar, NuZero, reynolds, Le, Bm, ambientTemperature, surfaceTemperature, I, writeData, vapor, calMethod)  

%% select the correlation
mikhil = strcmp(calMethod, 'mikhil');

%% select the fluid
waterVapor = strcmp(vapor, 'waterVapor');
ureaVapor = strcmp(vapor, 'ureaVapor');
nHeptaneVapor = strcmp(vapor, 'nHeptaneVapor');
nDecaneVapor = strcmp(vapor, 'nDecaneVapor');
vaporMixture = strcmp(vapor,'vaporMix');

%% Start of calculations

NuStar = zeros(2,1) ; 

if mikhil == 1
    
    CpGas = CpFilm ; 
    
    BtStart= 0.001;
    Bt = 0.1;

    if waterVapor == 1 || nHeptaneVapor == 1 || nDecaneVapor == 1 || vaporMixture == 1
    
        while ( abs((Bt-BtStart)/Bt) > 1.0e-11 )
    
            BtStart = Bt; 
    
            NuStar(2,1) = modifyNusselt(NuZero, BtStart, reynolds, calMethod) ;
    
            phi = (CpV/CpGas)*(ShStaar/NuStar(2,1))*(1/Le);
    
            Bt = (1+Bm)^phi - 1 ;
  
        end
    
    end
        
    if ureaVapor == 1
    
        if ge(surfaceTemperature,406.15)
    
            BtStart= 0.001;

            Bt = 0.1; 

            while ( abs((Bt-BtStart)/Bt) > 1.0e-11 )
        
                BtStart = Bt; 
    
                NuStar(2,1) = modifyNusselt(NuZero, BtStart, reynolds, calMethod) ;
    
                phi = (CpV/CpG)*(ShStaar/NuStar(2,1))*(1/Le);
    
                Bt = (1+Bm)^phi - 1 ; 

            end
    
        else
        
            NuStar(2,1) = 0; 
            
            Bt = 0 ; 
    
        end
    
    end

end

%% write data to files

savePath = './Properties/nonDimen/';
saveDataToFile = strcmp(writeData,'yes');
if saveDataToFile == 1
    
    if waterVapor == 1
        filename = [savePath,'BtWater','_',num2str(ambientTemperature),'K.txt'];
        
    elseif ureaVapor == 1
        filename = [savePath,'BtUrea','_',num2str(ambientTemperature),'K.txt'];
        
    elseif nHeptaneVapor == 1
        filename = [savePath,'BtHeptane','_',num2str(ambientTemperature),'K.txt'];
        
    elseif nDecaneVapor ==1 
        filename = [savePath,'BtDecane','_',num2str(ambientTemperature),'K.txt'];
        
    elseif vaporMixture ==1 
        filename = [savePath,'BtMixture','_',num2str(ambientTemperature),'K.txt'];
        
    end
        
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\n', 'i', 'Bt');
    end    
    fprintf(fileID, '%d\t %f\n', I, Bt);
    fclose(fileID);
        
    if waterVapor == 1
        filename = [savePath,'modNuWater','_',num2str(ambientTemperature),'K.txt'];
        
    elseif ureaVapor == 1
        filename = [savePath,'modNuUrea','_',num2str(ambientTemperature),'K.txt'];
        
    elseif nHeptaneVapor == 1
        filename = [savePath,'modNuHeptane','_',num2str(ambientTemperature),'K.txt'];
        
    elseif nDecaneVapor ==1 
        filename = [savePath,'modNuDecane','_',num2str(ambientTemperature),'K.txt'];
        
    elseif vaporMixture ==1 
        filename = [savePath,'modNuMixture','_',num2str(ambientTemperature),'K.txt'];
        
    end
        
    fileID = fopen(filename,'a');
    if I == 1
        fprintf(fileID, '%s\t %s\t %s\n', 'i', 'modNu');
    end    
    fprintf(fileID, '%d\t %d\t %f\n', I, NuStar(1,1), NuStar(2,1));
    fclose(fileID);

end
