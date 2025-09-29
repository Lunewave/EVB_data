close all; clear all; clc;
rospath = 'U:\Falcon_Project\20250811_MaranaDroneTest_32.45159N_111.21090W_Heading95E_RLCDir_-7\Square_2.447GHz\PA001_3811_LRMR_FALCONDroneTest_Square_2025-08-11-09-27-46';
path = 'U:\Falcon_Project\20250811_MaranaDroneTest_32.45159N_111.21090W_Heading95E_RLCDir_-7\Square_2.447GHz\';
csv = 'Aug-11th-2025-11-07AM-Flight-Airdata.csv';


%% ROTATOR LOCATION
clear v;
save_figs = 1;
ref_lat = 32.45159;       % North is positive
ref_lon = -111.21090;     % West is negative
ref_direction = 95;       % 0 is north, 90 is east, 180 is south. This is the direction that the 0 degree azimuth antenna is pointing.
noise_level_test = 45;

t_offset = 1.4;

multiple_freqs = 0;

%% Load Data
C_Time = load([rospath '\Camera_Time.mat']).Time(:, 2)+5998.75+64*.05;
C_dir = [rospath '\CameraData'];
C_Frames = dir(C_dir);

L_Time = load([rospath '\Lidar_Time.mat']).Time(:, 2)+5998.75+64*.05;
L_dir = [rospath '\LidarData'];
L_Frames = dir(L_dir);

R_Time = load([rospath '\radar_Time.mat']).Time(:, 2)+5998.75+64*.05;
R_dir = [rospath '\radar'];
R_Frames = dir(R_dir);


if multiple_freqs
    data_freq = [2.4218 2.4312] ; % Frequency of test data signal in GHz
    test_cache = fullfile(testpath, [num2str(data_freq(1)) '-' num2str(data_freq(2)) 'GHz_cached_test_data.mat']);
else
    data_freq = 2.447; %Frequency of test data signal in GHz
    test_cache = fullfile(path, [num2str(data_freq) 'GHz_cached_test_data.mat']);
end
load(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex', 'num_files', 'numgoodframes');


load([path 'MyTimestamps.mat'])
temp = load([path 'MyTimestamps.mat']);
fn = fieldnames(temp);
antenna_time = temp.(fn{1});
antenna_time = posixtime(antenna_time) + t_offset;
load([path 'AF_ITP_results_minh.mat'])
AF_ITP_results_minh(:, 1) = -AF_ITP_results_minh(:, 1);



% angle_offset = rad2deg(atan2(data.y(1), data.x(1)));
angle_offset = mod(90 - ref_direction, 360);
AF_ITP_results_minh(:, 1) = mod(AF_ITP_results_minh(:, 1) + angle_offset +180, 360) - 180;




%% Plot Figures
close all
l_pos = [0 0.8 -0.1];
r_pos = [0 0.8 0.06];
l_direction  = -7;
r_direction = -7;
LidarGT = nan(length(AF_ITP_results_minh(:, 1)), 4);

used_idx = []; % store lidar frame indices that have been used already

for ifile = 1:length(AF_ITP_results_minh(:, 1))

    % filter candidate indices so we don’t reuse old ones
    available_idx = setdiff(1:length(L_Time), used_idx);

    diff = antenna_time(ifile) - L_Time(available_idx);
    nonnegative_idx = find(diff >= 0);
    [a, Irel] = min(diff(nonnegative_idx));
    I = available_idx(nonnegative_idx(Irel));

    used_idx(end+1) = I; % mark this lidar frame as used

    L_Frame = load([L_dir '\Frame_' num2str(I) '.mat']);
    [xtemp, ytemp] = deal(cosd(l_direction)*L_Frame.X - sind(l_direction)*L_Frame.Y + l_pos(1), ...
                          sind(l_direction)*L_Frame.X + cosd(l_direction)*L_Frame.Y + l_pos(2));
    L_Frame = [xtemp, ytemp, L_Frame.Z + l_pos(3)];

    good_idx = find(L_Frame(:, 1) < 60);
    L_Frame = L_Frame(good_idx, :);
    good_idx = find(abs(L_Frame(:, 2)) < 30);
    L_Frame = L_Frame(good_idx, :);
    good_idx = find(L_Frame(:, 3) > -1);
    L_Frame = L_Frame(good_idx, :);


    % check if has points above 2m
    drone_det = find(L_Frame(:, 3) > 0.3 & L_Frame(:, 1) > 12 & L_Frame(:, 1) < 40 & ~(L_Frame(:, 1) >= 31.5 & L_Frame(:, 1) <= 32.5 & L_Frame(:, 2) >= -6.5   & L_Frame(:, 2) <= -4));
    % figure
    % subplot(1, 2, 1)
    % scatter3(L_Frame(:, 1), L_Frame(:, 2), L_Frame(:, 3))
    % hold on
    % plot(0:50, tand(AF_ITP_results_minh(ifile, 1)) * (0:50))
    % xlim([0 60])
    % ylim([-30 30])
    % zlim([-1 20])
    % view([0 90])
    % title('XY')
    % subplot(1, 2, 2)
    % scatter3(L_Frame(:, 1), L_Frame(:, 2), L_Frame(:, 3))
    % hold on
    % plot(0:50, tand(AF_ITP_results_minh(ifile, 1)) * (0:50))
    % view([0 0])
    % xlim([0 60])
    % ylim([-30 30])
    % zlim([-1 20])
    % title('XZ')
    % set(gcf, 'position', [100 , 100, 1400, 800])
    % sgtitle(['DF: ' num2str(ifile) ';   Lidar:' num2str(I)])

    while isempty(drone_det) && a < 0.25
        % remove this I from available set (without deleting from L_Time)
        used_idx(end+1) = I;

        available_idx = setdiff(1:length(L_Time), used_idx);
        [a, I] = min(abs(antenna_time(ifile) - L_Time(available_idx)));
        I = available_idx(I);

        used_idx(end+1) = I;

        L_Frame = load([L_dir '\Frame_' num2str(I) '.mat']);
        [xtemp, ytemp] = deal(cosd(l_direction)*L_Frame.X - sind(l_direction)*L_Frame.Y + l_pos(1), ...
                              sind(l_direction)*L_Frame.X + cosd(l_direction)*L_Frame.Y + l_pos(2));
        L_Frame = [xtemp, ytemp, L_Frame.Z + l_pos(3)];

        good_idx = find(L_Frame(:, 1) < 60);
        L_Frame = L_Frame(good_idx, :);
        good_idx = find(abs(L_Frame(:, 2)) < 30);
        L_Frame = L_Frame(good_idx, :);
        good_idx = find(L_Frame(:, 3) > -1);
        L_Frame = L_Frame(good_idx, :);

        drone_det = find(L_Frame(:, 3) > 0.3 & L_Frame(:, 1) > 12 & L_Frame(:, 1) < 40 & ~(L_Frame(:, 1) >= 31.5 & L_Frame(:, 1) <= 32.5 & L_Frame(:, 2) >= -6.5   & L_Frame(:, 2) <= -4));

        % figure
        % subplot(1, 2, 1)
        % scatter3(L_Frame(:, 1), L_Frame(:, 2), L_Frame(:, 3))
        % hold on
        % plot(0:50, tand(AF_ITP_results_minh(ifile, 1)) * (0:50))
        % xlim([0 60])
        % ylim([-30 30])
        % zlim([-1 20])
        % view([0 90])
        % title('XY')
        % subplot(1, 2, 2)
        % scatter3(L_Frame(:, 1), L_Frame(:, 2), L_Frame(:, 3))
        % hold on
        % plot(0:50, tand(AF_ITP_results_minh(ifile, 1)) * (0:50))
        % view([0 0])
        % xlim([0 60])
        % ylim([-30 30])
        % zlim([-1 20])
        % title('XZ')
        % set(gcf, 'position', [100 , 100, 1400, 800])
        % sgtitle(['DF: ' num2str(ifile) ';   Lidar:' num2str(I)])
    end
    if ~isempty(drone_det)
        disp('Detection')
        LidarGT(ifile, 4) = L_Time(I);
        LidarGT(ifile, 1:3) = L_Frame(drone_det(1), :);
    end
end



LidarGT = nan(length(L_Time), 4);

for ifile = 1:length(L_Time)


    L_Frame = load([L_dir '\Frame_' num2str(ifile) '.mat']);
    [xtemp, ytemp] = deal(cosd(l_direction)*L_Frame.X - sind(l_direction)*L_Frame.Y + l_pos(1), ...
                          sind(l_direction)*L_Frame.X + cosd(l_direction)*L_Frame.Y + l_pos(2));
    L_Frame = [xtemp, ytemp, L_Frame.Z + l_pos(3)];

    good_idx = find(L_Frame(:, 1) < 60);
    L_Frame = L_Frame(good_idx, :);
    good_idx = find(abs(L_Frame(:, 2)) < 30);
    L_Frame = L_Frame(good_idx, :);
    good_idx = find(L_Frame(:, 3) > -1);
    L_Frame = L_Frame(good_idx, :);


    % check if has points above 2m
    drone_det = find(L_Frame(:, 3) > 0.3 & L_Frame(:, 1) > 12 & L_Frame(:, 1) < 40 & ~(L_Frame(:, 1) >= 31.5 & L_Frame(:, 1) <= 32.5 & L_Frame(:, 2) >= -6.5   & L_Frame(:, 2) <= -4));

    if ~isempty(drone_det)
        disp('Detection')
        LidarGT(ifile, 4) = L_Time(ifile);
        LidarGT(ifile, 1:3) = L_Frame(drone_det(1), :);

    else
        disp('No Detection')
        LidarGT(ifile, 4) = L_Time(ifile);
        LidarGT(ifile, 1:3) = [NaN NaN NaN];
    end

end




%%
lidar_az = atan2d(LidarGT(:, 2), LidarGT(:, 1));
lidar_el = atan2d(LidarGT(:, 3), sqrt(LidarGT(:, 1).^2 + LidarGT(:, 2).^2));
% figure(101)
% tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% nexttile
% scatter(sqrt(LidarGT(:, 1).^2 + LidarGT(:, 2).^2 + LidarGT(:, 3).^2), abs(lidar_az - AF_ITP_results_minh(:, 1)));
% xlabel('Distance ($\sqrt{x^2 + y^2 + z^2}$) (m)', 'Interpreter', 'latex')
% ylabel('Azimuth Error abs(DF - Lidar GT) (deg)', 'Interpreter','latex')
% legend(['Mean: ' num2str(mean(abs(lidar_az - AF_ITP_results_minh(:, 1)), 'omitnan'))])
% nexttile
% scatter(sqrt(LidarGT(:, 1).^2 + LidarGT(:, 2).^2 + LidarGT(:, 3).^2), abs(lidar_el - AF_ITP_results_minh(:, 2)));
% xlabel('Distance ($\sqrt{x^2 + y^2 + z^2}$) (m)', 'Interpreter', 'latex')
% ylabel('Elevation Error abs(DF - Lidar GT) (deg)', 'Interpreter','latex')
% legend(['Mean: ' num2str(mean(abs(lidar_el - AF_ITP_results_minh(:, 2)), 'omitnan'))])
% set(gcf, 'position', [100 , 100, 1400, 800])
% sgtitle('Angle Error vs Distance', 'Interpreter', 'latex')


for i = 1:length(L_Time)

    [~, I] = min(abs(antenna_time - L_Time(i)));
    az = AF_ITP_results_minh(I, 1);

    if abs(lidar_az(i) - az) > 10
        lidar_az(i) = NaN;
        lidar_el(i) = NaN;
    end
end





figure(102)
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
plot(antenna_time - antenna_time(1), AF_ITP_results_minh(:, 1), '-o')
hold on
plot(LidarGT(:, 4) - antenna_time(1), lidar_az)
legend('DF Results', 'Lidar Results', 'Location', 'best')
xlabel('Time (s)', 'Interpreter', 'latex')
ylabel('Azimuth Angle (deg)', 'Interpreter', 'latex')
title('Azimuth Angle vs Time', 'Interpreter', 'latex')
nexttile
plot(antenna_time - antenna_time(1), AF_ITP_results_minh(:, 2), '-o')
hold on
plot(LidarGT(:, 4) - antenna_time(1), lidar_el)
legend('DF Results', 'Lidar Results', 'Location', 'best')
xlabel('Time (s)', 'Interpreter', 'latex')
ylabel('Elevation Angle (deg)', 'Interpreter', 'latex')
title('Elevation Angle vs Time', 'Interpreter', 'latex')
set(gcf, 'position', [100 , 100, 1400, 800])
sgtitle('Angle vs Time', 'Interpreter', 'latex')




%% Save Figures
if save_figs
    folderName = [num2str(data_freq*1000),'_MHz'];
    path = fullfile(path, folderName);
    
    if ~exist(path, 'dir')
        mkdir(path);
    else
        fprintf('Folder "%s" already exists.\n', folderName);
    end
    saveas(figure(101), fullfile(path, 'Error_vs_Lidar.jpeg'));
    saveas(figure(102), fullfile(path, 'DF_vs_Lidar_time.jpeg'));



end

