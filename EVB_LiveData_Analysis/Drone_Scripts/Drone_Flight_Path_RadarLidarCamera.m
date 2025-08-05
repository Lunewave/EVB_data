close all; clear all; clc;


%% ROTATOR LOCATION

save_figs = 0;
video = 0;
data_freq = 2.447; %Frequency of test data signal in GHz
ref_lat = 32.451553;       % North is positive
ref_lon = -111.211068;     % West is negative
ref_direction = 92;       % 0 is north, 90 is east, 180 is south. This is the direction that the 0 degree azimuth antenna is pointing.
noise_level_test = 45;




%% Load Data
rospath = uigetdir('U:\Falcon_Project\', 'Select Processed ROS bag data');
C_Time = load([rospath '\Camera_Time.mat']).Time(:, 2)+7*60*60-3.5;
C_dir = [rospath '\CameraData'];
C_Frames = dir(C_dir);

L_Time = load([rospath '\Lidar_Time.mat']).Time(:, 2)+7*60*60-3.5;
L_dir = [rospath '\LidarData'];
L_Frames = dir(L_dir);

R_Time = load([rospath '\radar_Time.mat']).Time(:, 2)+7*60*60-3.5;
R_dir = [rospath '\radar'];
R_Frames = dir(R_dir);


[csv, path] = uigetfile('U:\Falcon_Project\*.csv', 'Select CSV Flight Record');


test_cache = fullfile(path, [num2str(data_freq) 'GHz_cached_test_data.mat']);
load(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex', 'num_files', 'numgoodframes');


load([path 'MyTimestamps.mat'])
temp = load([path 'MyTimestamps.mat']);
fn = fieldnames(temp);
antenna_time = temp.(fn{1});
antenna_time = posixtime(antenna_time);
load([path 'AF_ITP_results.mat'])

fullpath = fullfile(path, csv);
T = readtable(fullpath);
data = struct();
for col = 1:width(T)
    data.(T.Properties.VariableNames{col}) = T{:, col};
end

data.datetime_utc_.TimeZone = 'UTC';
data.UTC_seconds = data.datetime_utc_ + seconds(mod(data.time_millisecond_, 1000)/1000);
data.UTC_seconds.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
data.UTC_seconds = posixtime(data.UTC_seconds);



ref_alt_feet = 1450;
ref_alt_m = ref_alt_feet * 0.3048; % Convert to meters

% Use your lat/lon and ref_lat/ref_lon and ref_alt_m
data.z = (data.altitude_above_seaLevel_feet_  - data.altitude_above_seaLevel_feet_(1))* 0.3048;
[data.x, data.y, ~] = geodetic2enu_custom(data.latitude, data.longitude, data.z, ref_lat, ref_lon, ref_alt_m);

good_idx = find(abs(data.x) < 1000);
data.x = data.x(good_idx);
data.y = data.y(good_idx);
data.z = data.z(good_idx);
data.UTC_seconds = data.UTC_seconds(good_idx);
% Smooth Time Data

% Assuming you have vectors x and y
p = polyfit(1:length(data.UTC_seconds), data.UTC_seconds, 1);         % Fit y = p(1)*x + p(2)
data.UTC_seconds = polyval(p, 1:length(data.UTC_seconds)) - 3.5;        % Evaluate the fitted line at x

dt = data.UTC_seconds - data.UTC_seconds(1);


% angle_offset = rad2deg(atan2(data.y(1), data.x(1)));
angle_offset = mod(90 - ref_direction, 360);
AF_ITP_results(:, 1) = mod(AF_ITP_results(:, 1) + angle_offset +180, 360) - 180;




%% Setup video writer
if video
    videoName = fullfile(path, 'Drone_Flight_Path_XY_XZ.mp4');
    v = VideoWriter(videoName, 'MPEG-4');
    v.FrameRate = 10;
    open(v);
    
    %% XY + XZ Video
    numPoints = length(data.UTC_seconds);
    trailLength = 5;  % last 5 points in trail
    baseSize = 36;    % largest marker size
    
    [sX, sY, sZ] = sphere;
    antenna_height = 2;
    l_pos = [0 -0.7 0];
    r_pos = [0 -1.43 0];
    l_direction  = -4;
    r_direction = -8;
    
    figure;
    set(gcf, 'WindowState', 'maximized');   % fill the screen completely
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Pre-compute limits
    margin = 10;
    xlimVals = [min([min(data.x), min(data.y)])-margin, max([max(data.x), max(data.y)])+margin];
    ylimVals = [min([min(data.x), min(data.y)])-margin, max([max(data.x), max(data.y)])+margin];
    zlimVals = [min([min(data.x), min(data.y)])-margin, max([max(data.x), max(data.y)])+margin];
    
    % Subplot 1: XY View
    axXY = nexttile;
    xlabel('X (m)'); ylabel('Y (m)');
    title('XY View');
    grid on; axis equal;
    xlim(xlimVals); ylim(ylimVals);
    hold on;
    surf(axXY, sX, sY, sZ + antenna_height);
    
    % Subplot 2: XZ View
    axXZ = nexttile;
    xlabel('X (m)'); ylabel('Z (m)');
    title('XZ View');
    grid on; axis equal;
    xlim(xlimVals); ylim(zlimVals);
    hold on;
    surf(axXZ, sX, sY + antenna_height, sZ);
    
    % Subplot 3: Radar/Lidar
    axRL = nexttile;
    xlabel('X (m)'); ylabel('Y (m)')
    title('Radar and Lidar')
    grid on; axis equal;
    xlim(xlimVals); ylim(ylimVals);
    hold on;
    
    % Subplot 4: Camera view
    axCAM = nexttile;
    axis off;
    title('Camera View');
    C_Frame = imread([C_dir '\Frame_1.jpg']);
    h_camera = imshow(C_Frame, 'Parent', axCAM);   % <- imshow used here ONCE
    
    
    
    % Initialize trails
    h_trailXY = gobjects(trailLength,1);
    h_trailXZ = gobjects(trailLength,1);
    for i = 1:trailLength
        h_trailXY(i) = plot(axXY, nan, nan, 'bo', 'MarkerFaceColor', 'b');
        h_trailXZ(i) = plot(axXZ, nan, nan, 'bo', 'MarkerFaceColor', 'b');
    end
    
    h_markerXY = plot(axXY, nan, nan, 'xk', 'LineWidth', 2, 'MarkerSize', 10);
    h_markerXZ = plot(axXZ, nan, nan, 'xk', 'LineWidth', 2, 'MarkerSize', 10);
    h_beamXY = plot(axXY, nan, nan, 'r--', 'LineWidth', 2);
    h_beamXZ = plot(axXZ, nan, nan, 'r--', 'LineWidth', 2);
    h_text = text(axXY, nan, nan, '', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);
    
    % LIDAR and RADAR point clouds (as scatter3 projections)
    h_lidarXY = plot(axRL, nan, nan, 'k.', 'MarkerSize', 1);
    % h_lidarXZ = plot(axRL, nan, nan, 'k.', 'MarkerSize', 2);
    h_radarXY = plot(axRL, nan, nan, 'm.', 'MarkerSize', 10);
    % h_radarXZ = plot(axRL, nan, nan, 'm.', 'MarkerSize', 10);
    
    
    for k = 1:numPoints
        % Find matching antenna frame
        [a, I] = min(abs(antenna_time - data.UTC_seconds(k)));
        if a<4
            curr_az = deg2rad(AF_ITP_results(I, 1));
            curr_el = deg2rad(AF_ITP_results(I, 2));
            r = sqrt(data.x(k)^2 + data.y(k)^2 + data.z(k)^2)*1.1;
            x = r*cos(curr_el)*cos(curr_az);
            y = r*cos(curr_el)*sin(curr_az);
            z = r*sin(curr_el) + antenna_height;
        else
            curr_az = NaN;
            curr_el = NaN;
            r = sqrt(data.x(k)^2 + data.y(k)^2 + data.z(k)^2)*1.1;
            x = r*cos(curr_el)*cos(curr_az);
            y = r*cos(curr_el)*sin(curr_az);
            z = r*sin(curr_el) + antenna_height;
        end
    
    
        %%% DO THIS FOR LIDAR, CAMERA, RADAR
    
        [~, I] = min(abs(L_Time - data.UTC_seconds(k)));
        L_Frame = load([L_dir '\Frame_' num2str(I) '.mat']);
        [xtemp ytemp] = deal(cosd(l_direction)*L_Frame.X - sind(l_direction)*L_Frame.Y + l_pos(1), sind(l_direction)*L_Frame.X + cosd(l_direction)*L_Frame.Y + l_pos(2));
        L_Frame = [xtemp, ytemp, L_Frame.Z+l_pos(3)];
    
    
        [~, I] = min(abs(R_Time - data.UTC_seconds(k)));
        R_Frame = load([R_dir '\Frame_' num2str(I) '.mat']);
        [xtemp ytemp] = deal(cosd(r_direction)*R_Frame.X - sind(r_direction)*R_Frame.Y + r_pos(1), sind(r_direction)*R_Frame.X + cosd(r_direction)*R_Frame.Y + r_pos(2));
        R_Frame = [xtemp, ytemp, R_Frame.Z+r_pos(3)];
    
    
        [~, I] = min(abs(C_Time - data.UTC_seconds(k)));
        C_Frame = imread([C_dir '\Frame_' num2str(I) '.jpg']);
    
    
        % --- Update LIDAR (green) ---
        set(h_lidarXY, 'XData', L_Frame(:,1), 'YData', L_Frame(:,2));
        % set(h_lidarXZ, 'XData', L_Frame(:,1), 'YData', L_Frame(:,3));
        
        % --- Update RADAR (magenta) ---
        set(h_radarXY, 'XData', R_Frame(:,1), 'YData', R_Frame(:,2));
        % set(h_radarXZ, 'XData', R_Frame(:,1), 'YData', R_Frame(:,3));
        
        % --- Update CAMERA view ---
        set(h_camera, 'CData', C_Frame);
    
    
    
    
        % Update beam
        set(h_beamXY, 'XData', [0, x], 'YData', [0, y]);
        set(h_beamXZ, 'XData', [0, x], 'YData', [antenna_height, z]);
    
        % Trail
        trailIdx = max(1, k-trailLength+1):k;
        for i = 1:trailLength
            if i <= length(trailIdx)
                idx = trailIdx(i);
                sizeFactor = baseSize / (length(trailIdx) - i + 1);
                set(h_trailXY(i), 'XData', data.x(idx), 'YData', data.y(idx), 'MarkerSize', sizeFactor/6);
                set(h_trailXZ(i), 'XData', data.x(idx), 'YData', data.z(idx), 'MarkerSize', sizeFactor/6);
            else
                set(h_trailXY(i), 'XData', nan, 'YData', nan);
                set(h_trailXZ(i), 'XData', nan, 'YData', nan);
            end
        end
    
        % Update drone marker
        set(h_markerXY, 'XData', data.x(k), 'YData', data.y(k));
        set(h_markerXZ, 'XData', data.x(k), 'YData', data.z(k));
    
        % Altitude label in XY plot
        set(h_text, 'Position', [data.x(k)+1, data.y(k)-1], ...
            'String', sprintf('Alt: %.1f m', data.z(k)));
    
    
        uistack(h_beamXY, 'top');
        uistack(h_beamXZ, 'top');
        uistack(h_markerXY, 'top');
        uistack(h_markerXZ, 'top');
        uistack(h_trailXY, 'top');
        uistack(h_trailXZ, 'top');
        % Range window around drone
        delta = 30;
        
        % Update XY plot limits
        xlim(axXY, [data.x(k) - delta, data.x(k) + delta]);
        ylim(axXY, [data.y(k) - delta, data.y(k) + delta]);
        
        % Update XZ plot limits
        xlim(axXZ, [data.x(k) - delta, data.x(k) + delta]);
        ylim(axXZ, [data.z(k) - delta, data.z(k) + delta]);
    
        % Update RL plot limits
        xlim(axRL, [data.x(k) - delta, data.x(k) + delta]);
        ylim(axRL, [data.z(k) - delta, data.z(k) + delta]);
    
        
        if k == 1
            legend(axXY, [h_trailXY(end), h_beamXY], {'Drone Position (GPS Estimation)', 'Antenna DF Results'}, 'Location', 'northeast');
            legend(axXZ, [h_trailXZ(end), h_beamXZ], {'Drone Altitude', 'Antenna DF Results'}, 'Location', 'northeast');
            legend(axRL, [h_lidarXY, h_radarXY], {'Lidar', 'Radar'}, 'Location', 'northeast');
        end
        
    
    
    
    
        drawnow;
        pause(0.1);
        
        frame = getframe(gcf);
        writeVideo(v, frame);
        sgtitle(['t = ' num2str(dt(k)) ' s'])
    end
    
    close(v);
    fprintf('Video saved to: %s\n', videoName);
end



%% Plot Figures
drone_az = rad2deg(atan2(data.y, data.x));
drone_el = rad2deg(atan2(data.z, sqrt(data.x.^2 + data.y.^2)));
drone_r = sqrt(data.x.^2 + data.y.^2 + data.z.^2);


for ifile = 1:length(AF_ITP_results(:, 1))
    [a, I] = min(abs(antenna_time(ifile) - data.UTC_seconds));
    if a<1 %Check to see if the drone time is close enough to the antenna time for good error.
        error(ifile, :) = [AF_ITP_results(ifile, 1) - drone_az(I), AF_ITP_results(ifile, 2) - drone_el(I)];
        drone_loc(ifile, :) = [drone_az(I) drone_el(I) drone_r(I)];
        drone_rxy(ifile) = sqrt(data.x(I).^2 + data.y(I).^2);

    else
        error(ifile, :) = [NaN NaN];
        drone_loc(ifile, :) = [NaN NaN NaN];
        drone_rxy(ifile) = NaN;

    end
end


figure(10)
subplot(1, 2, 1);
plot(error(:, 1), 'o-')
xlabel('File')
ylabel('Azimuth Error (degrees)')
title('Azimuth Error vs File')
subplot(1, 2, 2)
plot(error(:, 2), 'o-')
xlabel('File')
ylabel('Elevation Error (degrees)')
title('Elevation Error vs File')
set(gcf, 'Position', [100, 100, 1400, 700]);


figure(11)
subplot(1, 2, 1)
plot(data.UTC_seconds - data.UTC_seconds(1), drone_az, 'Color', [0, 0, 0]);
hold on
plot(antenna_time - data.UTC_seconds(1), AF_ITP_results(:, 1), 'o-')
xlabel('Seconds')
ylabel('Azimuth Angle (degrees)')
legend('Drone Ground Truth', 'FALCON DF Results', 'location', 'best')
subplot(1, 2, 2)
plot(data.UTC_seconds - data.UTC_seconds(1), drone_el, 'Color', [0, 0, 0]);
hold on
plot(antenna_time - data.UTC_seconds(1), AF_ITP_results(:, 2), 'o-')
xlabel('Seconds')
ylabel('Elevation Angle (degrees)')
legend('Drone Ground Truth', 'FALCON DF Results', 'location', 'best')
set(gcf, 'Position', [100, 100, 1400, 700]);


figure(12)
sgtitle(['6 Antenna SNR (Noise Level = ' num2str(noise_level_test) ' dB)'])
for i = 1:6
    subplot(2, 3, i)
    scatter(drone_loc(:, 3), Test_Mag(i, :)-noise_level_test)

    [~, sorted_idx] = sort(drone_loc(:, 3));
    drone_distance_sort = drone_loc(sorted_idx, 3);
    SNR = (Test_Mag(i, :)-noise_level_test)';
    SNR_sorted = SNR(sorted_idx);

    hold on
    X = 1 ./ drone_distance_sort;        % compute 1/r
    Y = SNR_sorted;
    validIdx = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
    X = X(validIdx);
    Y = Y(validIdx);


    p = polyfit(X, Y, 1);
    Y_fit = polyval(p, X);

    plot(1./X, Y_fit)
    eqnStr = sprintf('Fit: P = %.2f / r + %.3f', p(1), p(2));
    plot(NaN, NaN, 'w')
    ylabel('Power')
    xlabel('Distance (m)')
    ylim([0 2.5*10^4])
    title(['Antenna ' num2str(i)])
    grid on
    legend('Raw Data', '1/r Fit', eqnStr, 'Location', 'best')


    ylabel('SNR (dB)')
    xlabel('Distance (m)')
    ylim([0 55])
    title(['Antenna ' num2str(i)])
    grid on
end
set(gcf, 'Position', [100, 100, 1400, 700]);

figure(13)
sgtitle('6 Antenna Power Level')
for i = 1:6
    subplot(2, 3, i)
    scatter(drone_loc(:, 3), 10.^(Test_Mag(i, :)/20))

    [~, sorted_idx] = sort(drone_loc(:, 3));
    drone_distance_sort = drone_loc(sorted_idx, 3);
    mag_sorted = 10.^(Test_Mag(i, sorted_idx)/20)';

    hold on
    X = 1 ./ drone_distance_sort.^2;        % compute 1/r^2
    Y = mag_sorted;
    validIdx = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
    X = X(validIdx);
    Y = Y(validIdx);


    p = polyfit(X, Y, 1);
    Y_fit = polyval(p, X);

    plot(sqrt(1./X), Y_fit)
    eqnStr = sprintf('Fit: P = %.2f / r^{2} + %.3f', p(1), p(2));

    plot(NaN, NaN, 'w')
    ylabel('Power')
    xlabel('Distance (m)')
    ylim([0 2.5*10^4])
    title(['Antenna ' num2str(i)])
    grid on
    legend('Raw Data', '1/r^{2} Fit', eqnStr, 'Location', 'best')
end
set(gcf, 'Position', [100, 100, 1400, 700]);



figure(14)
sgtitle('Error vs Distance')
subplot(1, 2, 1)
scatter(drone_loc(:, 3), abs(error(:, 1)))
xlabel('Distance (m)')
ylabel('Absolute Azimuth Error (degrees)')
grid on
subplot(1, 2, 2)
scatter(drone_loc(:, 3), abs(error(:, 2)))
xlabel('Distance (m)')
ylabel('Absolute Elevation Error (degrees)')
grid on
set(gcf, 'Position', [100, 100, 1400, 700]);

figure(15)
sgtitle('Error vs Ground Truth')
subplot(1, 2, 1)
scatter(drone_loc(:, 1), error(:, 1))
xlabel('Ground Truth Drone Azimuth (degrees)')
ylabel('Azimuth Error (DF - GT) (degrees)')
xlim([-180 180])
grid on
subplot(1, 2, 2)
scatter(drone_loc(:, 2), error(:, 2))
xlabel('Ground Truth Drone Elevation (degrees)')
ylabel('Elevation Error (DF - GT) (degrees)')
xlim([0 90])
grid on
set(gcf, 'Position', [100, 100, 1400, 700]);


figure(16)
scatter(drone_rxy'.*cosd(AF_ITP_results(:, 1)), drone_rxy'.*sind(AF_ITP_results(:, 1)), 'filled')
hold on
for i = 1:num_files - 1
    x1 = drone_rxy(i)'*cosd(AF_ITP_results(i, 1));
    y1 = drone_rxy(i)'*sind(AF_ITP_results(i, 1));
    x2 = drone_rxy(i+1)'*cosd(AF_ITP_results(i+1, 1));
    y2 = drone_rxy(i+1)'*sind(AF_ITP_results(i+1, 1));
    line([x1 x2], [y1 y2], 'Color', [1 0 0 0.3], 'LineWidth', 1);  % Simulated transparency
    vec = [x2 - x1, y2 - y1];
    vec = vec / norm(vec);  % Normalize
    perp = [-vec(2), vec(1)];  % Perpendicular
    L = 0.4;  % Length of arrowhead
    W = 0.2;  % Width of arrowhead
    base = [x2, y2] - L * vec;
    arrow_x = [x2, base(1) + W*perp(1), base(1) - W*perp(1)];
    arrow_y = [y2, base(2) + W*perp(2), base(2) - W*perp(2)];
    patch(arrow_x, arrow_y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
for i = 1:num_files
    text(drone_rxy(i)'*cosd(AF_ITP_results(i, 1)), drone_rxy(i)'*sind(AF_ITP_results(i, 1)), num2str(i), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8)
end
xlabel('X (m)')
ylabel('Y (m)')
title('Signal Path')
grid on
set(gcf, 'Position', [100, 100, 1400, 700]);







%% Save Figures
if save_figs
    folderName = [num2str(data_freq*1000),'_MHz'];
    path = fullfile(path, folderName);
    
    if ~exist(path, 'dir')
        mkdir(path);
    else
        fprintf('Folder "%s" already exists.\n', folderName);
    end
    saveas(figure(10), fullfile(path, 'Error.jpeg'));
    saveas(figure(11), fullfile(path, 'Full_Flight.jpeg'));
    saveas(figure(12), fullfile(path, 'SNR.jpeg'));
    saveas(figure(13), fullfile(path, 'Power.jpeg'));
    saveas(figure(14), fullfile(path, 'Error_vs_Distance.jpeg'));
    saveas(figure(15), fullfile(path, 'Error_vs_GT.jpeg'));
    saveas(figure(16), fullfile(path, 'Path.jpeg'));


end

