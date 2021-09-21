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

%% This function calculates the drag coefficient.   
%  The equation may be suitably modified to take other expressions for the drag coefficient. 

function [dragCoefficient] = calDragCoefficient(Reynolds, ambientTemperature, I, writeData, calMethod)

saveDataToFile = strcmp(writeData,'yes');
savePath = './Properties/nonDimen/';

%% select the correlation
mikhil = strcmp(calMethod, 'mikhil');

%% calculate drag coefficient

if mikhil == 1
    
    dragCoefficient = (24/Reynolds) * (1 + (1/6) * Reynolds^(2/3)) ; 
    
    % write data to file
    if saveDataToFile == 1
        
        filename = [savePath,'Cd','_',num2str(ambientTemperature),'K.txt'];
              
        fileID = fopen(filename,'a');
        if I == 1
            fprintf(fileID, '%s\t %s\n', 'i', 'Cd');
        end    
        fprintf(fileID, '%d\t %f\n', I, dragCoefficient);
        fclose(fileID);
        
    end
    
end

