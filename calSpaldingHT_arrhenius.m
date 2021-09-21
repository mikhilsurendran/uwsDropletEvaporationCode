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

%% This function calculates the Spalding heat transfer number if an Arrhenius equation is used to estimate the depletion rate of urea 
%  Birkhold et al. (2007) UWS evaporation modeling

function [Bt, phi, NuStar] = calSpaldingHT_arrhenius(CpV, CpFilm, vaporizationEnthalpy, ShStaar, NuZero, reynolds, Le, ambientTemperature, surfaceTemperature, boilingPointOne, I, writeData, vapor, calMethod)  

%% select the correlation
mikhil = strcmp(calMethod, 'mikhil');

%% select the fluid
waterVapor = strcmp(vapor, 'waterVapor');
ureaVapor = strcmp(vapor, 'ureaVapor');
nHeptaneVapor = strcmp(vapor, 'nHeptaneVapor');
nDecaneVapor = strcmp(vapor, 'nDecaneVapor');
vaporMixture = strcmp(vapor,'vaporMix');

%% Start of calculations

if mikhil == 1
     
    Bt = CpFilm * (ambientTemperature - surfaceTemperature)/vaporizationEnthalpy ;
    
    Bt_BoilPt = CpFilm * (ambientTemperature - boilingPointOne)/vaporizationEnthalpy ;
    
    NuStar(1,1) = modifyNusselt(NuZero, Bt_BoilPt, reynolds, calMethod) ; % used to calculate the evaporation rate of water after it reaches the boiling point
    NuStar(2,1) = modifyNusselt(NuZero, Bt, reynolds, calMethod) ;        % use to calculate heat transfer rate 
    
    phi = (CpV/CpFilm)*(ShStaar/NuStar(2,1))*(1/Le) ;

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
        fprintf(fileID, '%s\t %s\t %s\n', 'i', 'modNu', 'modNu');
    end    
    fprintf(fileID, '%d\t %d\t %f\n', I, NuStar(1,1), NuStar(2,1));
    fclose(fileID);

end