function [CL,CD,CM,wind_aoa,k, wind_rho,no_wind_N_force] = BoeingClean(wing_area,BoeingChord)
k = 0;
for i = [9 10 11 12 21 22 23 24 33 34 35 36]
    k = k + 1;
% All Boeing Cases (Only studying the clean case)
if(i >= 10)
    Boeing = strcat('data/787_G',num2str(i),'.csv');
else
    Boeing = strcat('data/787_G0',num2str(i),'.csv');
end

% Load File
file = load(Boeing);

% Store data from all no wind cases
no_wind_pressure = file(1:300,1); % Find the atmosph. pressure
no_wind_A_force = file(1:300,25); % Find the axial force acting on the sting
no_wind_N_force = file(1:300,24); % Find normal force acting on the sting
no_wind_aoa = file(1:300,23); % Find angle of attack from file
no_wind_pitch_moment = file(1:300,26); % Find pitching moment acting on the sting
no_wind_v_inf = file(1:300,4);
no_wind_rho = file(1:300,3);
distance = 0.063; % Distance from the sting in m

% Store data from all wind cases
wind_pressure(k,:) = file(301:600,1); % Find the atmosph. pressure
wind_A_force = file(301:600,25); % Find the axial force acting on the sting
wind_N_force = file(301:600,24); % Find the normal force acting on the sting
wind_aoa(k,:) = file(301:600,23); % Find the angle of attack
wind_pitch_moment = file(301:600,26); % Find the pitching moment acting on the sting
wind_v_inf(k,:) = file(301:600,4);
wind_rho = file(301:600,3);
wind_q_inf(k,:) = file(301:600,5);

% Calculate Aerodynamic Properties from Data
Normal(k,:) = wind_N_force - no_wind_N_force; % Normal Force
mat = Normal(k,:)*distance; % Moment about the sting to subtract
Axial(k,:) = wind_A_force - no_wind_A_force; % Axial Force
Moment(k,:) = wind_pitch_moment - no_wind_pitch_moment - transpose(mat); % Pitching Moment about COG

Lift(k,:) = Normal(k,:).*cosd(wind_aoa(k,:))-Axial(k,:).*sind(wind_aoa(k,:)); % Lift Force
Drag(k,:) = Normal(k,:).*sind(wind_aoa(k,:))+Axial(k,:).*cosd(wind_aoa(k,:)); % Drag Force
%Alpha = linspace(-8, 20, 299);

% Calculate CL, CD, and CM
CL(k, :) = Lift(k,:) ./ (wind_q_inf(k, :) * wing_area);
CD(k, :) = Drag(k,:) ./ (wind_q_inf(k, :) * wing_area);
CM(k,:) = Moment(k,:) ./ (wind_q_inf(k, :) * wing_area * BoeingChord);

end


end