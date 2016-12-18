function [CL,CD,CM,wind_aoa,k,wind_rho,no_wind_N_force,std_q,error_CL,error_CD] = Cleanload(wing_area,F_16_chord)
%% All "Clean" Cases
k = 0;
for i = [1 2 3 4 13 14 15 16 25 26 28]
    k = k + 1;
if(i >= 10)
    f_16 = strcat('data/F16_Clean_G',num2str(i),'.csv');
else
    f_16 = strcat('data/F16_Clean_G0',num2str(i),'.csv');
end

%% Load Data File
file = load(f_16);

%% Store data from all no wind cases
no_wind_pressure = file(1:300,1); % Find the atmosph. pressure
no_wind_A_force = file(1:300,25); % Find the axial force acting on the sting
no_wind_N_force = file(1:300,24); % Find normal force acting on the sting
no_wind_aoa = file(1:300,23); % Find angle of attack from file
no_wind_pitch_moment = file(1:300,26); % Find pitching moment acting on the sting
no_wind_v_inf = file(1:300,4);
no_wind_rho = file(1:300,3);
distance = 0.0144; % Distance from the sting in m

%% Store data from all wind cases
wind_pressure(k,:) = file(301:600,1); % Find the atmosph. pressure
wind_A_force = file(301:600,25); % Find the axial force acting on the sting
wind_N_force = file(301:600,24); % Find the normal force acting on the sting
wind_aoa(k,:) = file(301:600,23); % Find the angle of attack
wind_pitch_moment = file(301:600,26); % Find the pitching moment acting on the sting
wind_v_inf(k,:) = file(301:600,4); % Wind speed
wind_rho(k,:) = file(301:600,3); % Wind air density
wind_q_inf(k,:) = file(301:600,5); % Wind Dynamic Pressure 

%% Calculate Aerodynamic Properties from Data
Normal(k,:) = wind_N_force - no_wind_N_force; % Normal Force
mat = Normal(k,:)*distance;
Axial(k,:) = wind_A_force - no_wind_A_force; % Axial Force
Moment(k,:) = wind_pitch_moment- no_wind_pitch_moment - transpose(mat); % Pitching Moment

Lift(k,:) = Normal(k,:).*cosd(wind_aoa(k,:))-Axial(k,:).*sind(wind_aoa(k,:)); % Lift Force
Drag(k,:) = Normal(k,:).*sind(wind_aoa(k,:))+Axial(k,:).*cosd(wind_aoa(k,:)); % Drag Force
%Alpha = linspace(-8, 20, 299);

%% Calculate CL, CD, and CM
CL(k, :) = Lift(k,:) ./ (wind_q_inf(k, :) * wing_area);
CD(k, :) = Drag(k,:) ./ (wind_q_inf(k, :) * wing_area);
CM(k,:) = Moment(k,:) ./ (wind_q_inf(k, :) * wing_area * F_16_chord);

% Calculate Standard Deviation of data for comparison
std_alpha_clean = std(wind_aoa(k,:));
std_normal_clean = std(Normal(k,:));
std_axial_clean = std(Axial(k,:));
std_moment_clean = std(Moment(k,:));
std_v_inf = mean(std(wind_v_inf(k,:)));

std_q(k, :) = sqrt((wind_rho(k,:).*wind_v_inf(k, :)).^2.*std_v_inf.^2);

std_lift(k, :) = sqrt(cosd(wind_aoa(k, :)).^2.*std_normal_clean.^2 + ...
        (sind(wind_aoa(k, :)).^2).*std_axial_clean.^2 + ...
        (Normal(k, :).*sind(wind_aoa(k, :)) - ...
        Axial(k, :).*cosd(wind_aoa(k, :))).^2.*std_alpha_clean);
    
    std_drag(k, :) = sqrt(sind(wind_aoa(k, :)).^2.*std_normal_clean.^2 + ...
        (cosd(wind_aoa(k, :)).^2).*std_axial_clean.^2 + ...
        (Normal(k, :).*cosd(wind_aoa(k, :)) + ...
        Axial(k, :).*sind(wind_aoa(k, :))).^2.*std_alpha_clean);
    
    error_CL(k, :) = sqrt(((1./(wind_q_inf(k, :).* ...
        wing_area)).^2).*std_lift(k, :).^2 + ...
        (-Lift(k,:)./(wind_q_inf(k, :).*wing_area).^2).^2.*std_q(k, :).^2);   
    
    error_CD(k, :) = sqrt(((1./(wind_q_inf(k, :).* ...
        wing_area)).^2).*std_drag(k, :).^2 + ...
        (-Drag(k,:)./(wind_q_inf(k, :).*wing_area).^2).^2.*std_q(k, :).^2);
end

end