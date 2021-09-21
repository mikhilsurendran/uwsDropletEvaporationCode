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

%% This function plots all data. 

function plotData(time, diameter, initialDiameter, MFtwo, temperature, squareOfDia, squareOfDiaNormalized, tNormalized, totalMass, massOfComponentOne, massOfComponentTwo, changeInMass, changeInMassOne, changeInMassTwo, boilPt, isItUWS)

%% estimate the array size

[numRows, ~] = size(time); 

%% Find the time at which the temperature of the droplet becomes equal to the boiling point of water. 

if ge(max(temperature), boilPt)
    
    for count = 1:numRows
        
        if ge(temperature(count,1), boilPt)
            
            timeAtBoilPt = time(count,1) ;
            
            break ;
            
        end
        
    end
    
end    

%% Find the time at which the temperature of the droplet becomes equal to the melting point of urea. 

if strcmp(isItUWS, 'urea') == 1 
    
    if ge(max(temperature), 406.15)
    
        for count = 1:numRows
        
            if ge(temperature(count,1), 406.15)
            
                timeAtMeltPt = time(count,1) ;
            
                break ;
            
            end
        
        end
    
    end
    
else 
    
    timeAtMeltPt = 0 ;
    
end

%% Find the time at which mass fraction of the second component becomes 100 % 

for count = 1:numRows
    
    if MFtwo(count,1) > 0.99999 
            
        timeAtUnity = time(count,1) ;
            
        break ;
            
    end
    
end


%% start plotting

figure(1)
yyaxis left;

plot (time, diameter, 'b','LineWidth',2);
ylim ([0 ((initialDiameter*1000)*1.2)]);
ylabel('Droplet Diameter (mm)');

hold on; 
yyaxis right;

plot (time, MFtwo, 'black','LineWidth',2);
grid on ; 
grid minor ;
ylim ([0 1.2]);
ylabel('Component Two Mass Fraction');

y= ylim ; 
if ge(max(temperature), boilPt)   % plot a vertical line when droplet temperature = boiling point
    hold on;    
    plot([timeAtBoilPt timeAtBoilPt],[y(1) y(2)],'--m') ;
    text(timeAtBoilPt, 0.1, 'T = BP', 'HorizontalAlignment', 'right');
end
if ge(max(temperature), 406.15)   % plot a vertical line when droplet temperature = melting point of urea 
    hold on;      
    plot([timeAtMeltPt timeAtMeltPt],[y(1) y(2)],'--m') ;
    text(timeAtMeltPt, 0.2, 'T = MP', 'HorizontalAlignment', 'right');
end
if max(MFtwo) > 0.99999    % plot a vertical line when component two mass fraction = 1
    hold on;
    plot([timeAtUnity timeAtUnity],[y(1) y(2)],'--m') ;
    text(timeAtUnity, 0.1, 'Y = 1', 'HorizontalAlignment', 'left');
end

legend('diameter', 'component two mass fraction', 'Location', 'northeast');

xlabel('Time (s)');

% ----------------------------------------------------------------------------------------------------------------------------------------------------

figure(2)
yyaxis left;

plot (time, temperature, 'r','LineWidth',2);
ylim ([200 (max(temperature * 1.2))]);
ylabel('Droplet Temperature (K)');

hold on; 
yyaxis right;

plot (time, MFtwo, 'black','LineWidth',2);
grid on ; 
grid minor ;
ylim ([0 1.2]);
ylabel('Component Two Mass Fraction / (D/Do)^2');

hold on; 
yyaxis right;

plot (time, squareOfDiaNormalized, 'blue','LineWidth',2);

y= ylim ; 
if ge(max(temperature), boilPt)   % plot a vertical line when droplet temperature = boiling point
    hold on; 
    y= ylim ;     
    plot([timeAtBoilPt timeAtBoilPt],[y(1) y(2)],'--m') ;
    text(timeAtBoilPt, 0.1, 'T = BP', 'HorizontalAlignment', 'right');
end
if ge(max(temperature), 406.15)   % plot a vertical line when droplet temperature = melting point of urea 
    hold on;      
    plot([timeAtMeltPt timeAtMeltPt],[y(1) y(2)],'--m') ;
    text(timeAtMeltPt, 0.2, 'T = MP', 'HorizontalAlignment', 'right');
end
if max(MFtwo) > 0.99999   % plot a vertical line when component two mass fraction = 1
    hold on; 
    plot([timeAtUnity timeAtUnity],[y(1) y(2)],'--m') ;
    text(timeAtUnity, 0.1, 'Y = 1', 'HorizontalAlignment', 'left');
end

legend('temperature', 'component two mass fraction', '(d/do)^2', 'Location', 'northeast');

xlabel('Time (s)');

% ----------------------------------------------------------------------------------------------------------------------------------------------------

figure(3)
yyaxis left;

plot (time, squareOfDia, 'b','LineWidth',2);
ylim ([0 (((initialDiameter*1000)^2)*1.2)]);
ylabel('Droplet Diameter^2 (mm^2)');

hold on; 
yyaxis right;

plot (time, MFtwo, 'black','LineWidth',2);
ylim ([0 1.2]);
ylabel('Component Two Mass Fraction');

y= ylim ; 
if ge(max(temperature), boilPt)   % plot a vertical line when droplet temperature = boiling point
    hold on; 
    y= ylim ;     
    plot([timeAtBoilPt timeAtBoilPt],[y(1) y(2)],'--m') ;
    text(timeAtBoilPt, 0.1, 'T = BP', 'HorizontalAlignment', 'right');
end
if ge(max(temperature), 406.15)   % plot a vertical line when droplet temperature = melting point of urea 
    hold on;      
    plot([timeAtMeltPt timeAtMeltPt],[y(1) y(2)],'--m') ;
    text(timeAtMeltPt, 0.2, 'T = MP', 'HorizontalAlignment', 'right');
end
if max(MFtwo) > 0.99999   % plot a vertical line when component two mass fraction = 1 
    hold on;
    plot([timeAtUnity timeAtUnity],[y(1) y(2)],'--m') ;
    text(timeAtUnity, 0.1, 'Y = 1', 'HorizontalAlignment', 'left');
end

legend('d2', 'component two mass fraction', 'Location', 'northeast');

xlabel('Time (s)');

% ----------------------------------------------------------------------------------------------------------------------------------------------------

figure(4)
yyaxis left;

plot (tNormalized, squareOfDiaNormalized, 'b','LineWidth',2); 
ylim ([0 1.2]);
ylabel('(D/Do)^2');

hold on; 
yyaxis right;

plot (tNormalized, MFtwo, 'black','LineWidth',2);
grid on ; 
grid minor ;
ylim ([0 1.2]);
ylabel('Component Two Mass Fraction');

legend('(d/do)^2', 'component two mass fraction', 'Location', 'northeast');

xlabel('Normalized Time (s/mm^2)');

% ----------------------------------------------------------------------------------------------------------------------------------------------------

figure(5)
yyaxis left;

plot (time, squareOfDiaNormalized, 'b','LineWidth',2); 
ylim ([0 1.2]);
ylabel('(D/Do)^2');

hold on; 
yyaxis right;

plot (time, MFtwo, 'black','LineWidth',2);
grid on ; 
grid minor ;
ylim ([0 1.2]);
ylabel('Component Two Mass Fraction');

y= ylim ; 
if ge(max(temperature), boilPt)   % plot a vertical line when droplet temperature = boiling point
    hold on; 
    y= ylim ;     
    plot([timeAtBoilPt timeAtBoilPt],[y(1) y(2)],'--m') ;
    text(timeAtBoilPt, 0.1, 'T = BP', 'HorizontalAlignment', 'right');
end
if ge(max(temperature), 406.15)   % plot a vertical line when droplet temperature = melting point of urea 
    hold on;      
    plot([timeAtMeltPt timeAtMeltPt],[y(1) y(2)],'--m') ;
    text(timeAtMeltPt, 0.2, 'T = MP', 'HorizontalAlignment', 'right');
end
if max(MFtwo) > 0.99999    % plot a vertical line when component two mass fraction = 1
    hold on;
    plot([timeAtUnity timeAtUnity],[y(1) y(2)],'--m') ;
    text(timeAtUnity, 0.1, 'Y = 1', 'HorizontalAlignment', 'left');
end

legend('(d/do)^2', 'component two mass fraction', 'Location', 'northeast');

xlabel('Time (s)');

% ----------------------------------------------------------------------------------------------------------------------------------------------------

figure(6)
yyaxis left;

plot (time, totalMass, 'b','LineWidth',2); 
ylabel('mass (ug)');

hold on; 

plot (time, massOfComponentOne, '--g','LineWidth',2); 

hold on; 

plot (time, massOfComponentTwo, '--r','LineWidth',2); 


hold on; 
yyaxis right;

plot (time, MFtwo, 'black','LineWidth',2);
grid on ; 
grid minor ;
ylim ([0 1.2]);
ylabel('Component Two Mass Fraction');

y= ylim ; 
if ge(max(temperature), boilPt)   % plot a vertical line when droplet temperature = boiling point
    hold on; 
    y= ylim ;     
    plot([timeAtBoilPt timeAtBoilPt],[y(1) y(2)],'--m') ;
    text(timeAtBoilPt, 0.1, 'T = BP', 'HorizontalAlignment', 'right');
end
if ge(max(temperature), 406.15)   % plot a vertical line when droplet temperature = melting point of urea 
    hold on;      
    plot([timeAtMeltPt timeAtMeltPt],[y(1) y(2)],'--m') ;
    text(timeAtMeltPt, 0.2, 'T = MP', 'HorizontalAlignment', 'right');
end
if max(MFtwo) > 0.99999    % plot a vertical line when component two mass fraction = 1
    hold on;
    plot([timeAtUnity timeAtUnity],[y(1) y(2)],'--m') ;
    text(timeAtUnity, 0.1, 'Y = 1', 'HorizontalAlignment', 'left');
end

legend('total mass', 'component one mass', 'component two mass', 'component two mass fraction', 'Location', 'northeast');

xlabel('Time (s)');

% ----------------------------------------------------------------------------------------------------------------------------------------------------

figure(7)
yyaxis left;

plot (time, totalMass, 'b','LineWidth',2); 
ylabel('mass (ug)');

hold on; 

plot (time, massOfComponentOne, '--g','LineWidth',2);

hold on; 

plot (time, massOfComponentTwo, '--r','LineWidth',2); 



hold on; 
yyaxis right;

plot (time, temperature, 'black','LineWidth',2);
grid on ; 
grid minor ;
ylim ([200 (max(temperature * 1.2))]);
ylabel('Droplet Temperature (K)');

y= ylim ; 
if ge(max(temperature), boilPt)   % plot a vertical line when droplet temperature = boiling point 
    hold on; 
    y= ylim ;     
    plot([timeAtBoilPt timeAtBoilPt],[y(1) y(2)],'--m') ;
    text(timeAtBoilPt, 225, 'T = BP', 'HorizontalAlignment', 'right');
end
if ge(max(temperature), 406.15)   % plot a vertical line when droplet temperature = melting point of urea 
    hold on;      
    plot([timeAtMeltPt timeAtMeltPt],[y(1) y(2)],'--m') ;
    text(timeAtMeltPt, 215, 'T = MP', 'HorizontalAlignment', 'right');
end
if max(MFtwo) > 0.99999    % plot a vertical line when component two mass fraction = 1
    hold on; 
    plot([timeAtUnity timeAtUnity],[y(1) y(2)],'--m') ;
    text(timeAtUnity, 225, 'Y = 1', 'HorizontalAlignment', 'left');
end

legend('total mass', 'component one mass', 'component two mass', 'temperature', 'Location', 'northeast');

xlabel('Time (s)');

% ----------------------------------------------------------------------------------------------------------------------------------------------------

figure(8)
yyaxis left;

plot (time, changeInMass, 'b','LineWidth',2); 
ylabel('rate of change of mass (ug/s)');

hold on; 

plot (time, changeInMassOne, '--g','LineWidth',2); 

hold on; 

plot (time, changeInMassTwo, '--r','LineWidth',2); 


hold on; 
yyaxis right;

plot (time, MFtwo, 'black','LineWidth',2);
grid on ; 
grid minor ;
ylim ([0 1.2]);
ylabel('Component Two Mass Fraction');

y= ylim ; 
if ge(max(temperature), boilPt)   % plot a vertical line when droplet temperature = boiling point
    hold on; 
    y= ylim ;     
    plot([timeAtBoilPt timeAtBoilPt],[y(1) y(2)],'--m') ;
    text(timeAtBoilPt, 0.1, 'T = BP', 'HorizontalAlignment', 'right');
end
if ge(max(temperature), 406.15)   % plot a vertical line when droplet temperature = melting point of urea 
    hold on;      
    plot([timeAtMeltPt timeAtMeltPt],[y(1) y(2)],'--m') ;
    text(timeAtMeltPt, 0.2, 'T = MP', 'HorizontalAlignment', 'right');
end
if max(MFtwo) > 0.99999    % plot a vertical line when component two mass fraction = 1
    hold on;
    plot([timeAtUnity timeAtUnity],[y(1) y(2)],'--m') ;
    text(timeAtUnity, 0.1, 'Y = 1', 'HorizontalAlignment', 'left');
end

legend('change in droplet mass', 'change in component one mass', 'change in component two mass', 'component two mass fraction', 'Location', 'northeast');

xlabel('Time (s)');

% ----------------------------------------------------------------------------------------------------------------------------------------------------

figure(9)
yyaxis left;

plot (time, changeInMass, 'b','LineWidth',2); 
ylabel('rate of change of mass (ug/s)');

hold on; 

plot (time, changeInMassOne, '--g','LineWidth',2); 

hold on; 

plot (time, changeInMassTwo, '--r','LineWidth',2); 


hold on; 
yyaxis right;

plot (time, temperature, 'black','LineWidth',2);
grid on ; 
grid minor ;
ylim ([200 (max(temperature * 1.2))]);
ylabel('Droplet Temperature (K)');

y= ylim ; 
if ge(max(temperature), boilPt)   % plot a vertical line when droplet temperature = boiling point
    hold on; 
    y= ylim ;     
    plot([timeAtBoilPt timeAtBoilPt],[y(1) y(2)],'--m') ;
    text(timeAtBoilPt, 225, 'T = BP', 'HorizontalAlignment', 'right');
end
if ge(max(temperature), 406.15)   % plot a vertical line when droplet temperature = melting point of urea 
    hold on;      
    plot([timeAtMeltPt timeAtMeltPt],[y(1) y(2)],'--m') ;
    text(timeAtMeltPt, 215, 'T = MP', 'HorizontalAlignment', 'right');
end
if max(MFtwo) > 0.99999     % plot a vertical line when component two mass fraction = 1
    hold on; 
    plot([timeAtUnity timeAtUnity],[y(1) y(2)],'--m') ;
    text(timeAtUnity, 225, 'Y = 1', 'HorizontalAlignment', 'left');
end

legend('change in droplet mass', 'change in component one mass', 'change in component two mass', 'temperature', 'Location', 'northeast');

xlabel('Time (s)');

