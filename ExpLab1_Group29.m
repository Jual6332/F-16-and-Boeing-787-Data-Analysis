%% Original Author: Justin Alvey
% Lab: Experimental Lab 2004 - F-16 and 787 Windtunnel Exp.
% Date Created: 1/14/16 
% Date Last Edited: 3/1/16

clc;
clear all;
close all;

%% Test Data 
test1 = 'F16_CLEAN_LA.csv';
test2 = 'F16_CLEAN_LA2.csv';

%% Constants for Calculation
scaleF16 = 1/48; % Scale of the F16
scale787 = 1/225; % Scale of the 787
rho_Boulder = 1.0425; % Density at Boulder (kg/m^3)
pound2N = 4.44822162; % Conversion from lb to N

%% F16 Wing Area
wing_area_actual = 27.87; % Reference: http://www.lockheedmartin.com/us/products/f16/F-16Specifications.html
wing_area = wing_area_actual*(scaleF16^2); % Wing area sized-down
WeightF16 = 133446.64800000002; % Landing weight of the full-scale F16 (N)

%% F16 Chord
% Cited From "Modern Combat Aircraft Design" by Klaus Huenecke: 
F_16_root_chord = 5.04; % Actual length
F_16_tip_chord =  1.0583; % Actual length
F_16_chord = scaleF16*(F_16_root_chord + F_16_tip_chord)/2; % Scaled-down chord (F16)

%% Boeing 787 Chord
% Cited From: http://www.lissys.demon.co.uk/samp1/
root_787_chord = 38.94; % Actual Length
tip_787_chord = 5.55; % Actual Length
BoeingChord = scale787*(root_787_chord + tip_787_chord)/2; % Scaled-down chord (787)

%% 787 Wing Area
% Cited From: https://en.wikipedia.org/wiki/General_Dynamics_F-16_Fighting_Falcon
wing_area_787_actual = 325; % Area of the 787 wing
wing_area_787 = wing_area_787_actual*(scale787^2); % Scaled-down wing area of the 787
Weight787 = 1601359.776; % Landing weight of the Boeing 787 (N)

%% Load F16 Clean Data Files from Directory (Current Folder)
[CLclean, CDclean, CMclean, wind_aoa_clean,k,wind_rho,weightmodelF16,error_CLclean, error_CDclean] = Cleanload(wing_area,F_16_chord); % Calls Cleanload.m

%% Load F16 Loaded Data Files from Directory (Current Folder)
[CLdirty, CDdirty, CMdirty, wind_aoa_dirty,c,error_CLdirty, error_CDdirty] = Dirtyload(wing_area,F_16_chord); % Calls Dirtyload.m

%% Load Boeing 787 Clean Data Files from Directory (Current Folder)
[CL787clean, CD787clean, CM787clean, wind_aoa_787clean,d,wind_rho_787,weightmodel787] = BoeingClean(wing_area_787,BoeingChord); 
% Calls BoeingClean.m

%% Error Calculations 
%Only use -5 to 5 degrees because past that the data is not as linear, so
%it will throw off the best fit line a little
%aoa = linspace(-5,5,11); % Angle of attack for polyfit calculations
%alpha = linspace(-8,20,29); % Angle of attack for average graphs

alpha_clean_std = std(wind_aoa_clean);
alpha_dirty_std = std(wind_aoa_dirty);

CM_clean_std = std(CMclean);
CM_dirty_std = std(CMclean);

CL_clean_std = std(CLclean);
CD_clean_std = std(CDclean);

CL_dirty_std = std(CLclean);
CD_dirty_std = std(CDdirty);

CL_787 = std(CL787clean);
CD_787 = std(CD787clean);
CM_787 = std(CM787clean);
alpha_787 = std(wind_aoa_787clean);

% Normalizing the Data for Each Standard Deviation
for i=0:14
    % Boeing 787
    CL_787_stds(i*20+1:20+i*20) = mean(CL_787(i*20+1:20+i*20));
    CD_787_stds(i*20+1:20+i*20) = mean(CD_787(i*20+1:20+i*20));
    CM_787_stds(i*20+1:20+i*20) = mean(CM_787(i*20+1:20+i*20));
    alpha_787_stds(i*20+1:20+i*20) = mean(alpha_787(i*20+1:20+i*20));
    
    alpha_clean_stds(i*20+1:20+i*20) = mean(alpha_clean_std(i*20+1:20+i*20));
    alpha_dirty_stds(i*20+1:20+i*20) = mean(alpha_dirty_std(i*20+1:20+i*20));
    
    % Coefficients of Lift
    CL_clean_stds(i*20+1:20+i*20) = mean(CL_clean_std(i*20+1:20+i*20));
    CLclstd(i+1) = CL_clean_stds(i*20+1);
    CL_dirty_stds(i*20+1:20+i*20) = mean(CL_dirty_std(i*20+1:20+i*20));
    
    % Coefficients of Drag
    CD_clean_stds(i*20+1:20+i*20) = mean(CD_clean_std(i*20+1:20+i*20));
    CD_dirty_stds(i*20+1:20+i*20) = mean(CD_dirty_std(i*20+1:20+i*20));
    
    % Coefficients of Moments
    CM_clean_stds(i*20+1:20+i*20) = mean(CM_clean_std(i*20+1:20+i*20));
    CM_dirty_stds(i*20+1:20+i*20) = mean(CM_dirty_std(i*20+1:20+i*20));
    
end

for i= [2:20 22:40 42:60 62:80 82:100 102:120 122:140 142:160 162:180 182:200 202:222 224:240 242:260 262:280 282:300]
    for j=1:15
    if CL_clean_stds(i) == CLclstd(j)
        % Boeing
        alpha_787_stds(i) = 0; 
        CL_787_stds(i) = 0;
        CD_787_stds(i) = 0;
        CM_787_stds(i) = 0;
        
        alpha_clean_stds(i) = 0;
        alpha_dirty_stds(i) = 0;
        CL_clean_stds(i) = 0;
        CL_dirty_stds(i) = 0;
        CD_clean_stds(i) = 0;
        CD_dirty_stds(i) = 0;
        CM_clean_stds(i) = 0;
        CM_dirty_stds(i) = 0;
        
    end
    end
end

% Locate the array elements in the std array that are not repeated
% Find returns the index of the values greater than 0, which are the singly
% repeated values

% Boeing Data
alpha787find = find(alpha_787_stds>0); 
CL787find = find(CL_787_stds>0);
CD787find = find(CD_787_stds>0);
CM787find = find(CM_787_stds>0);

% F-16 Data
alphacleanfind = find(alpha_clean_stds>0);
alphadirtyfind = find(alpha_dirty_stds>0);
CLcleanfind = find(CL_clean_stds>0);
CLdirtyfind =  find(CL_dirty_stds>0);
CDcleanfind =  find(CD_clean_stds>0);
CDdirtyfind =  find(CD_dirty_stds>0);
CMcleanfind = find(CM_clean_stds>0);
CMdirtyfind = find(CM_clean_stds>0);

% Reassign the non-repeated std values to a new array
finalalpha787std = alpha_787_stds(alpha787find);
finalCL787std = CL_787_stds(CL787find);
finalCD787std = CD_787_stds(CD787find);
finalCM787std = CM_787_stds(CM787find);

finalalphacleanstd = alpha_clean_stds(alphacleanfind);
finalalphadirtystd = alpha_dirty_stds(alphadirtyfind);
finalCLcleanstd = CL_clean_stds(CLcleanfind);
finalCLdirtystd = CL_dirty_stds(CLdirtyfind);
finalCDcleanstd = CD_clean_stds(CDcleanfind);
finalCDdirtystd = CD_dirty_stds(CDdirtyfind);
finalCMcleanstd = CM_clean_stds(CMcleanfind);
finalCMdirtystd = CM_dirty_stds(CMdirtyfind);


%% Aerodynamic Calculations for the Clean F16
fprintf('1. Stability Data for the Clean F16\n');

% Calculating Static Longitudinal Stability - Clean F16
% Assumption Made: the The relationship for CL vs alpha is "nearly" linear
% for the angles of attack used
coeff_clean = polyfit(wind_aoa_clean, CMclean,1); % Coefficients for the polyfit of the CM data
dCM = coeff_clean(1); % dCm/d_alpha when alpha = 0 is the first coefficient of the polyfit coefficients
fprintf('dCm/d(alpha) for Clean F16: %0.4f \n',dCM) % Format the values for the user

% Calculating dCL/d(alpha) for Clean F16
% Assumption Made: the The relationship for CM vs alpha is "nearly" linear
% for the angles of attack used
coeff_clean = polyfit(wind_aoa_clean,CLclean,1);
dCL = coeff_clean(1);
fprintf('dCL/d(alpha) for Clean F16: %0.4f \n',dCL) % Format the values for the user

% Calculating the Static Margin (-dCM/d(alpha) / dCL/d(alpha)) - Clean F16
staticMargin = -dCM/dCL;
fprintf('-(dCM/d(alpha))/(dCL/d(alpha)) for Clean F16: %0.4f \n',staticMargin) % Format the values for the user
fprintf('\n')

% Calculating the Minimum Landing Speed - Clean F16
Vstall = sqrt((2*max(weightmodelF16))./(mean(mean(wind_rho))*wing_area*max(max(CLclean))));
fprintf('2. Minimum Landing Speed of the Clean F16 Model: Vmin = %0.3f m/s \n',1.2*Vstall)
fprintf('\n')

%% Aerodynamic Calculations for the Dirty F16
fprintf('3. Stability Data for the Dirty F16\n');

% Calculating Static Longitudinal Stability - Dirty F16
coeff_clean = polyfit(wind_aoa_dirty, CMdirty,1); % Coefficients for the polyfit of the CM data
dCM = coeff_clean(1); % dCm/d_alpha when alpha = 0 is the first coefficient of the polyfit coefficients
fprintf('dCm/d(alpha) for Dirty F16: %0.4f \n',dCM) % Format the values for the user

% Calculating dCL/d(alpha) - Dirty F16
coeff_clean = polyfit(wind_aoa_dirty,CLdirty,1);
dCL = coeff_clean(1);
fprintf('dCL/d(alpha) for Dirty F16: %0.4f \n',dCL) % Format the values for the user

% Calculating the Static Margin (-dCM/d(alpha) / dCL/d(alpha)) - Dirty F16
staticMargin = -dCM/dCL;
fprintf('-(dCM/d(alpha))/(dCL/d(alpha)) for Dirty F16: %0.4f \n',staticMargin) % Format the values for the user
fprintf('\n')

% Calculating the Minimum Landing Speed - Dirty F16
Vstall = sqrt((2*max(weightmodelF16))./(mean(mean(wind_rho))*wing_area*max(max(CLdirty))));
fprintf('4. Minimum Landing Speed of the Dirty F16 Model: Vmin = %0.3f m/s \n',1.2*Vstall)
fprintf('\n')

% Calculating the Minimum Landing Speed of the Full-Scale F16
Vstall_fullsize = sqrt((2*WeightF16)./(rho_Boulder*wing_area_actual*max(max(CLclean))));
fprintf('5. Minimum Landing Speed of the Full-size F16: Vmin = %0.3f m/s \n',1.2*Vstall_fullsize)
fprintf('\n')

%% Aerodynamic Calculations for the Boeing 787
fprintf('6. Stability Data for the Boeing 787 \n');

% Calculating Static Longitudinal Stability - Boeing 787
coeff_dirty = polyfit(wind_aoa_787clean, CM787clean,1); % Coefficients for the polyfit of the CM data
dCM = coeff_dirty(1); % dCm/d_alpha when alpha = 0 is the first coefficient of the polyfit coefficients
fprintf('dCm/d(alpha) for Boeing 787: %0.4f \n',dCM) % Format the values for the user

% Calculating dCL/d(alpha) - Boeing 787
coeff_clean = polyfit(wind_aoa_787clean, CL787clean,1);
dCL = coeff_clean(1);
fprintf('dCL/d(alpha) for Boeing 787: %0.4f \n',dCL) % Format the values for the user

% Calculating the Static Margin (-dCM/d(alpha) / dCL/d(alpha)) - Boeing 787
staticMargin = -dCM/dCL;
fprintf('-(dCM/d(alpha))/(dCL/d(alpha)) for Boeing 787: %0.4f \n',staticMargin) % Format the values for the user
fprintf('\n')

% Calculating the Minimum Landing Speed of the Boeing 787
Vstall = sqrt((2*max(weightmodel787))./(mean(mean(wind_rho_787))*wing_area_787*max(max(CL787clean))));
fprintf('7. Minimum Landing Speed of the Model Boeing 787: Vmin = %0.3f m/s \n',1.3*Vstall)
fprintf('\n')

% Calculating the Minimum Landing Speed of the Actual Boeing 787
Vstall = sqrt((2*Weight787)./(rho_Boulder*wing_area_787_actual*max(max(CL787clean))));
fprintf('8. Minimum Landing Speed of the Actual Boeing 787: Vmin = %0.3f m/s \n',1.3*Vstall)
fprintf('\n')

% Calculate the Lift/Drag ratio of the F-16
CLF16avg = max(max(CLclean));
CDF16avg = max(max(CDclean));
lift_drag = CLF16avg/CDF16avg;% Calculate the lift-to-drag ratio 
fprintf('9. Lift/Drag Ratio of the F-16: L/D = %0.3f \n',lift_drag)
fprintf('\n')

% Calculate the Lift/Drag ratio of the Boeing 787
CL787avg = max(max(CL787clean));
CD787avg = max(max(CD787clean));
lift_drag = CL787avg/CD787avg;% Calculate the lift-to-drag ratio 
fprintf('10. Lift/Drag Ratio of the Boeing 787: L/D = %0.3f \n',lift_drag)
fprintf('\n')

%% Plotting - Clean F16
% Plot CL versus alpha before error analysis
figure(1);
   for j = 1:k
       hold on
       plot(wind_aoa_clean(j,:),CLclean(j,:));
       xlabel('Angle of attack (degrees)')
       ylabel('Coefficient of Lift')
       title('C_{L} vs \alpha for Clean F-16')
       xlim([min(min(wind_aoa_clean)) 21])
       ylim([-0.75 1.45])
   end

% Plot CD versus alpha before error analysis
figure(2);
   for j = 1:k
       hold on;
       plot(CDclean(j,:),wind_aoa_clean(j,:));
       xlabel('Angle of attack (degrees)')
       ylabel('Coefficient of Drag')
       title('C_{D} vs \alpha for Clean F-16')
       ylim([-10 21])
   end
   
   % Plot CD versus CL before error analysis
   figure(8);
   for j = 1:k
       hold on;
       plot(CDclean(j,:),CLclean(j,:));
       xlabel('Coefficient of Drag')
       ylabel('Cofficient of Lift')
       title('C_{D} vs C_{L} for Clean F-16')
       xlim([0 0.4])
       ylim([-0.75 1.45])
   end
   
   figure(3);   
   % Plot CM versus alpha before error analysis
   for j = 1:k
       hold on;
       plot(wind_aoa_clean(j,:),CMclean(j,:));
       xlabel('Angle ot attack (deg)')
       ylabel('Cofficient of Moment')
       title('C_{M} vs \alpha for Clean F-16')
       xlim([min(min(wind_aoa_clean)) 21])
       ylim([-0.15 0.135])
   end
   
%% Plotting - Dirty F16
% Plot CL versus alpha before error analysis
figure(4);
   for j = 1:c
       hold on
       plot(wind_aoa_dirty(j,:),CLdirty(j,:));
       xlabel('Angle of attack (degrees)')
       ylabel('Coefficient of Lift')
       title('C_{L} vs \alpha for "Dirty" F-16')
       xlim([min(min(wind_aoa_dirty)) 21])
       ylim([-0.8 1.5])
   end
   
figure(5);
% Plot CD versus alpha before error analysis
   for j = 1:c
       hold on;
       plot(CDdirty(j,:),wind_aoa_dirty(j,:));
       xlabel('Angle of attack (degrees)')
       ylabel('Coefficient of Drag')
       title('C_{D} vs \alpha for "Dirty" F-16')
       ylim([-10 21])
       xlim([0.03 0.6])
   end
   
   figure(6);
   % Plot CD versus CL before error analysis
   for j = 1:c
       hold on;
       plot(CDdirty(j,:),CLdirty(j,:));
       xlabel('Coefficient of Lift')
       ylabel('Cofficient of Drag')
       title('C_{D} vs C_{L} for "Dirty" F-16')
       xlim([0.03 0.6])
       ylim([-0.85 1.5])
   end
   
   figure(7);   
   % Plot CM versus alpha before error analysis
   for j = 1:c
       hold on;
       plot(wind_aoa_dirty(j,:),CMdirty(j,:));
       xlabel('Angle ot attack (deg)')
       ylabel('Cofficient of Moment')
       title('C_{M} vs \alpha for "Dirty" F-16')
       xlim([min(min(wind_aoa_dirty)) 21])
       ylim([-0.27 0.3])
   end
   
%% Plotting - Boeing 787
% Plot CL versus alpha before error analysis
figure(9);
for i=1:d
    hold on;
    plot(wind_aoa_787clean(i,:),CL787clean(i,:))
    xlabel('Angle of attack (deg)')
    ylabel('Coefficient of Lift')
    title('C_{L} vs \alpha for Boeing 787')
    xlim([min(min(wind_aoa_787clean)) 21])
    ylim([-0.75 1.5])
    
end

% Plot CD versus alpha before error analysis
figure(10);
for i=1:d
    hold on;
    plot(CD787clean(i,:),wind_aoa_787clean(i,:))
    xlabel('Angle of attack (deg)')
    ylabel('Coefficient of Drag')
    title('C_{D} vs \alpha for Boeing 787')
    ylim([-10,21])
    xlim([0 0.56])
    
end


% Plot CD versus CL before error analysis
figure(11);
for i=1:d
    hold on;
    plot(CL787clean(i,:),CD787clean(i,:))
    xlabel('Coefficient of Lift')
    ylabel('Coefficient of Drag')
    title('C_{D} vs C_{L} for Boeing 787')
    xlim([min(min(CL787clean))-0.05,1.5])
    ylim([0 0.56])
    
end

% Plot CM versus alpha before error analysis
figure(12);
for i=1:d
    hold on;
    plot(wind_aoa_787clean(i,:),CM787clean(i,:))
    xlabel('Angle of attack (deg)')
    ylabel('Coefficient of Moment')
    title('C_{M} vs \alpha for Boeing 787')
    %xlim([min(min(wind_aoa_787clean)),21])
    %ylim([-10 22])    
end

% CM vs alpha for both clean and dirty F-16s
figure(13);
for i=1:c
    hold on;
    plot(wind_aoa_clean(i,:),CMclean(i,:))
    plot(wind_aoa_dirty(i,:),CMdirty(i,:))
    xlabel('Angle of attack (deg)')
    ylabel('Coefficient of Moment')
    title('C_{M} vs \alpha for both Clean and Dirty F-16')
    xlim([min(min(wind_aoa_clean)),20.5])
end

%% Error Calculations
m = mean(CLclean);
n = mean(CLdirty);
o = mean(wind_aoa_clean);
p = mean(wind_aoa_dirty);

% Average CL vs alpha for Clean and Dirty Cases (2 stds)
figure
for i=1:c
hold on
xlabel('Angle of Attack (degrees)')
ylabel('Coefficient of Lift (C_L)')
title(['Average C_{L} vs \alpha with Two' ...
    ' Std - Unloaded F16'])
plot(p(CLdirtyfind),n(CLdirtyfind),'b','LineWidth',2)
errorbar(p(CLdirtyfind),n(CLdirtyfind), 2*finalCLcleanstd,'r')
herrorbar(p(CLdirtyfind),n(CLdirtyfind), 2*finalalphacleanstd,'r')
legend('Unloaded F16', '2 \sigma', 'Location', 'Northwest')
ylim([-0.75 1.45])
xlim([-10 22])
end

% Average CL vs alpha for Dirty Case (2 stds)
figure
for i=1:c
hold on
xlabel('Angle of Attack (degrees)')
ylabel('Coefficient of Lift (C_L)')
title(['Average C_{L} vs \alpha with Two' ...
    ' Std - Loaded F16'])
plot(p(CLdirtyfind),m(CLdirtyfind),'b','LineWidth',2)
errorbar(p(CLdirtyfind),m(CLdirtyfind), 2*finalCLdirtystd,'r')
herrorbar(p(CLdirtyfind),m(CLdirtyfind), 2*finalalphadirtystd,'r')
legend('Unloaded F16', '2 \sigma', 'Location', 'Northwest')
ylim([-0.75 1.45])
xlim([-10 22])
end

u = mean(CDclean);
v = mean(CDdirty);

% Average CD vs CL for Clean Case (2 stds)
figure
for i=1:c
hold on
xlabel('Coefficient of Lift (C_{L})')
ylabel('Coefficient of Drag (C_{D})')
title(['Average C_{D} vs Average C_{L} with Two' ...
    ' Std - Unloaded F16'])
plot(m(CDcleanfind),u(CDcleanfind),'b','LineWidth',2)
%plot(n(CDdirtyfind),v(CDdirtyfind),'b','LineWidth',2)
errorbar(m(CDcleanfind),u(CDcleanfind), 2*finalCDcleanstd,'r')
%errorbar(n(CDdirtyfind),v(CDdirtyfind), 2*finalCDdirtystd,'b')
herrorbar(m(CDcleanfind),u(CDcleanfind), 2*finalCLcleanstd,'r')
%herrorbar(n(CDdirtyfind),v(CDdirtyfind), 2*finalCLdirtystd,'b')
legend('Unloaded F16', '2 \sigma', 'Location', 'Northwest')
xlim([-0.75 1.5])
ylim([0 0.43])
end

% Average CD vs CL for Dirty Case (2 stds)
figure
for i=1:c
hold on
xlabel('Coefficient of Lift (C_{L})')
ylabel('Coefficient of Drag (C_{D})')
title(['Average C_{D} vs Average C_{L} with Two' ...
    ' Std - Loaded F16'])
plot(n(CDdirtyfind),v(CDdirtyfind),'b','LineWidth',2)
errorbar(n(CDdirtyfind),v(CDdirtyfind), 2*finalCDdirtystd,'r')
herrorbar(n(CDdirtyfind),v(CDdirtyfind), 2*finalCLdirtystd,'r')
legend('Loaded F16', '2 \sigma', 'Location', 'Northwest')
xlim([-0.75 1.5])
ylim([0.05 0.57])
end

f = mean(CMclean);
g = mean(CMdirty);

% Average CM vs alpha for Clean Case (2 stds)
figure
for i=1:c
hold on
xlabel('Angle of Attack (degrees)')
ylabel('Coefficient of Drag (C_{M})')
title(['Average C_{M} vs Average \alpha with Two' ...
    ' Std - Unloaded F16'])
plot(o(CMcleanfind),f(CMcleanfind),'b','LineWidth',2)
errorbar(o(CMcleanfind),f(CMcleanfind), 2*finalCMcleanstd,'r')
herrorbar(o(CMcleanfind),f(CMcleanfind), 2*finalalphacleanstd,'r')
legend('Unloaded F16', '2 \sigma', 'Location', 'Northwest')
ylim([-0.15 0.14])
xlim([-10 22])
end

% Average CM vs alpha for Dirty Case (2 stds)
figure
for i=1:c
hold on
xlabel('Angle of Attack (degrees)')
ylabel('Coefficient of Drag (C_{M})')
title(['Average C_{M} vs Average \alpha with Two' ...
    ' Std - Loaded F16'])
plot(p(CMcleanfind),g(CMcleanfind),'b','LineWidth',2)
errorbar(p(CMcleanfind),g(CMcleanfind), 2*finalCMdirtystd,'r')
herrorbar(p(CMcleanfind),g(CMcleanfind), 2*finalalphadirtystd,'r')
legend('Loaded F16', '2 \sigma', 'Location', 'Northwest')
ylim([-0.22 0.3])
xlim([-10 22])
end

%% Error Calculations - Boeing 787 

m = mean(wind_aoa_787clean);
o = mean(CL787clean);
p = mean(CD787clean);
q = mean(CM787clean);

% Average CL vs alpha for Boeing 787 (2 stds)
figure
for i=1:d
hold on
xlabel('Angle of Attack (degrees)')
ylabel('Coefficient of Lift (C_{L})')
title(['Average C_{L} vs Average \alpha with Two' ...
    ' Std - Boeing 787'])
plot(m(CL787find),o(CL787find),'b','LineWidth',2)
errorbar(m(CL787find),o(CL787find), finalCL787std,'r')
herrorbar(m(CL787find),o(CL787find), 2*finalalpha787std,'r')
legend('Boeing 787', '2 \sigma', 'Location', 'Northwest')
%ylim([-0.22 0.3])
xlim([-10 22])
end

% Average CD vs alpha for Boeing 787 (2 stds)
figure
for i=1:d
hold on
xlabel('Angle of Attack (degrees)')
ylabel('Coefficient of Drag (C_{D})')
title(['Average C_{D} vs Average C_{L} with Two' ...
    ' Std - Boeing 787'])
plot(o(CD787find),p(CD787find),'b','LineWidth',2)
errorbar(o(CD787find),p(CD787find), finalCD787std/4,'r')
herrorbar(o(CD787find),p(CD787find), finalalpha787std/4,'r')
legend('Boeing 787', '2 \sigma', 'Location', 'Northwest')
%ylim([-0.22 0.3])
%xlim([-10 22])
end

% Average CM vs alpha for Boeing 787 (2 stds)
figure
for i=1:d
hold on
xlabel('Angle of Attack (degrees)')
ylabel('Coefficient of Moment (C_{M})')
title(['Average C_{M} vs Average \alpha with Two' ...
    ' Std - Boeing 787'])
plot(m(CM787find),q(CM787find),'b','LineWidth',2)
errorbar(m(CM787find),q(CM787find), finalCM787std,'r')
herrorbar(m(CM787find),q(CM787find), 2*finalalpha787std,'r')
legend('Boeing 787', '2 \sigma', 'Location', 'Northeast')
%ylim([-0.22 0.3])
xlim([-10 22])
end

%% END OF CODE