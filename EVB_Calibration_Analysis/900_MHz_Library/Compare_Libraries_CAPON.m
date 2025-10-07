close all; clear all; clc;


AZ_start = 180; AZ_end = -180; AZ_step = -3;
EL_start = 66; EL_end = -6; EL_step = -3;
save_figs = 1;
lib_location = 'Calibration Library 915 MHz';
test_location = 'Test Calibration Library 915 MHz';
frequency = 915;
noise_level_cal = 45;

shifts1 = [0 0 0 0 0 0];
shifts2 = [0 0 0 0 0 0];

shifts11 = [0 0 0 0 0 0];
shifts12 = [0 0 0 0 0 0];


libpath = 'U:\Direction_Finding\20250930_MaranaCalibrationLibrary_915MHz_360AZ_66_to_-6EL';
testpath = 'U:\Direction_Finding\20250930_MaranaCalibrationLibrary_915MHz_360AZ_66_to_-6EL';

%%%%%%%%%%% FIXED PARAMETERS %%%%%%%%%%%%%
AZ_data = AZ_start:AZ_step:AZ_end;
AZ_steps = length(AZ_data);
EL_data = EL_start:EL_step:EL_end;
EL_steps = length(EL_data);
numpeaks2check = 1; %# of peaks to check in each dimension of angle interpolation
files = dir('C:\Users\WorkStation2\Documents\GitHub\EVB_data\EVB_Calibration_Analysis\900_MHz_Library\MUSIC\Result_0p915_lib_CAPON\*.mat');
AF_ITP_results = zeros(length(files), 2);
for i = 1:length(files)
    load(fullfile(files(i).folder, [num2str(i) '_0p915_result.mat']));
    AF_ITP_results(i, 1) = final_AZ_ITP;
    AF_ITP_results(i, 2) = final_EL_ITP;
end

AZ_results = AF_ITP_results(:, 1);
EL_results = AF_ITP_results(:, 2);


AZ_err_ITP = reshape(AZ_results, [121, 25]);
EL_err_ITP = reshape(EL_results, [121, 25]);


for i = 1:121
    EL_err_ITP(i, :) = EL_err_ITP(i, :) - (66:-3:-6);
end

for i = 1:25
    AZ_err_ITP(:, i) = AZ_err_ITP(:, i) - (180:-3:-180)';
end


AZ_err_ITP=mod(AZ_err_ITP+180,360)-180;
EL_err_ITP=mod(EL_err_ITP+180,360)-180;
%% Plot Angle Errors
%%%%%%%%%%%%%%%%%%%  angle error  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp=AZ_err_ITP;
AZ_err_ave_ITP=mean(abs(tmp(:)));
AZ_err_std_ITP=std(abs(tmp(:)));
AZ_err_max_ITP=max(abs(tmp(:)));

tmp=EL_err_ITP;
EL_err_ave_ITP=mean(abs(tmp(:)));
EL_err_std_ITP=std(abs(tmp(:)));
EL_err_max_ITP=max(abs(tmp(:)));

step_error = 10;

figure(10)
subplot(1, 2, 1)

imagesc(AZ_data,EL_data,abs(AZ_err_ITP.'));
caxis(abs(AZ_step) * [-0, step_error]);
cb = colorbar;
set(cb, 'Ticks', (-step_error * 0) : abs(AZ_step) : (step_error * abs(AZ_step)));
xlabel('AZ');ylabel('EL');

title(['Interpolated AZ error at ' num2str(frequency)  ' MHz']);
err_text = sprintf(['Average Interpolated Azimuth Error: %.3f degrees\n' ...
                    'Average Interpolated Azimuth Error SD: %.3f degrees\n' ...
                    'Max Interpolated Azimuth Error: %.3f degrees'], AZ_err_ave_ITP, AZ_err_std_ITP, AZ_err_max_ITP);
annotation('textbox', [0.01, 0.025, 0.3, 0.05], ...
    'String', err_text, ...
    'FitBoxToText', 'on', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 9);
subplot(1, 2, 2)

imagesc(AZ_data,EL_data,abs(EL_err_ITP.'));
caxis(abs(EL_step) * [-0, step_error]);
cb = colorbar;
set(cb, 'Ticks', (step_error * 0) : abs(EL_step) : (-step_error * EL_step));
xlabel('AZ');ylabel('EL');

title(['Interpolated EL error at ' num2str(frequency)  ' MHz']);
err_text = sprintf(['Average Interpolated Elevation Error: %.3f degrees\n' ...
                    'Average Interpolated Elevation Error SD: %.3f degrees\n' ...
                    'Max Interpolated Elevation Error: %.3f degrees'], EL_err_ave_ITP, EL_err_std_ITP, EL_err_max_ITP);
annotation('textbox', [0.47, 0.025, 0.3, 0.05], ...
    'String', err_text, ...
    'FitBoxToText', 'on', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 9);
set(gcf, 'Position', [100, 100, 1400, 700]);




figure(12)
subplot(1,2,1)
hold on;
plot(AZ_data,AZ_data.','o-');
plot(AZ_data,AZ_err_ITP(:,23).'+AZ_data,'o-');
xlabel('AZ (deg)');ylabel('AZ angle (deg)');
grid on; ylim([AZ_end-5 AZ_start+5]);xlim([AZ_end-5 AZ_start+5])
legend('Real Angle', 'Interpolated Angle', 'Location', 'best')
title('Azimuth at EL = 0')
subplot(1,2,2)
hold on;
plot([EL_start:EL_step:EL_end],[EL_start:EL_step:EL_end].','o-');
plot([EL_start:EL_step:EL_end],EL_err_ITP(61,:)+[EL_start:EL_step:EL_end],'o-');
xlabel('EL (deg)');ylabel('EL angle (deg)');
grid on; ylim([EL_end-5 EL_start+5]);xlim([EL_end-5 EL_start+5])
legend('Real Angle', 'Interpolated Angle', 'Location', 'best')
title('Elevation at AZ = 0')
set(gcf, 'Position',  [200, 200, 1200, 500]);

figure(13)
subplot(1, 2, 1)
plot(AZ_data,AZ_err_ITP(:,23).');
grid on
xlabel('Azimuth')
ylabel('Error')
title('Elevation  = 0 Azimuth Error')
subplot(1, 2, 2)
plot(EL_data,EL_err_ITP(61,:).');
grid on
xlabel('Elevation')
ylabel('Error')
title('Azimuth  = 0 Elevation Error')
set(gcf, 'Position',  [200, 200, 1200, 500]);






%% Save Figures
if save_figs
    folderName = [num2str(frequency),'_MHz'];
    newFolderPath = fullfile(testpath, folderName);
    
    if ~exist(newFolderPath, 'dir')
        mkdir(newFolderPath);
    else
        fprintf('Folder "%s" already exists.\n', folderName);
    end
    
    

    saveas(figure(10), fullfile(newFolderPath, 'Error_CAPON.jpeg'))

    saveas(figure(12), fullfile(newFolderPath, 'Principle_Plane_CAPON.jpeg'))
    
    saveas(figure(13), fullfile(newFolderPath, 'Principle_Plane_Error_CAPON.jpeg'))
end


clear a b c d e1
a = reshape(EL_err_ITP ,[1 121*25]);
b = reshape(AZ_err_ITP ,[1 121*25]);
figure()
plot(a)
grid on
xlabel('File #')
ylabel('EL Error')

c = a>6;
d = find(c == 1);
e1 = mod(d-1, 121)+1;
e2 = floor((d - 1) / 121) + 1;   % elevation index (increments every 121)
az_ang_axis = 180:-3:-180;
el_ang_axis = 66:-3:-6;

az_angs = az_ang_axis(e1);
el_angs = el_ang_axis(e2);
% Histogram with 1-wide bins
figure()
% Bin edges centered on measurement angles
bin_width = 3;  % since your AZ step is 3 degrees
edges = (min(az_ang_axis)-0.5) : 1 : (max(az_ang_axis)+0.5);
histogram(flip(az_angs), edges);
xlabel('AZ Angle');
ylabel('Count');
title('Histogram of AZ angle with EL Err > 6');




figure()
% Bin edges centered on measurement angles
bin_width = 3;  % since your AZ step is 3 degrees
edges = (min(el_ang_axis)-0.5) : 1 : (max(el_ang_axis)+0.5);
histogram(flip(el_angs), edges);
xlabel('EL Angle');
ylabel('Count');
title('Histogram of EL angle with EL Err > 6');




figure()
plot(b)
grid on
xlabel('File #')
ylabel('AZ Error')

c = abs(b)>50;
d = find(c == 1);
e1 = mod(d-1, 121)+1;
e2 = floor((d - 1) / 121) + 1;   % elevation index (increments every 121)
az_ang_axis = 180:-3:-180;
el_ang_axis = 66:-3:-6;

az_angs = az_ang_axis(e1);
el_angs = el_ang_axis(e2);
% Histogram with 1-wide bins
figure()
% Bin edges centered on measurement angles
bin_width = 3;  % since your AZ step is 3 degrees
edges = (min(az_ang_axis)-0.5) : 1 : (max(az_ang_axis)+0.5);
histogram(flip(az_angs), edges);
xlabel('AZ Angle');
ylabel('Count');
title('Histogram of AZ angle with AZ abs(Err) > 50');




figure()
% Bin edges centered on measurement angles
bin_width = 3;  % since your AZ step is 3 degrees
edges = (min(el_ang_axis)-0.5) : 1 : (max(el_ang_axis)+0.5);
histogram(flip(el_angs), edges);
xlabel('EL Angle');
ylabel('Count');
title('Histogram of EL angle with AZ abs(Err) > 50');