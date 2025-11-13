close all; clear all; clc;


%% ROTATOR LOCATION

save_figs = 1;
data_freq = 2.405; %Frequency of test data signal in GHz
ref_lat = 32.451886;       % North is positive
ref_lon = -111.210883;     % West is negative
ref_direction = 130;       % 0 is north, 90 is east, 180 is south. This is the direction that the 0 degree azimuth antenna is pointing.
noise_level_test = 45;

t_offset = 10.75;

%% Load Data
[csv, path] = uigetfile('U:\Direction_Finding\*.csv', 'Select CSV Flight Record');


test_cache = fullfile(path, [num2str(data_freq) 'GHz_cached_test_data.mat']);
load(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex', 'num_files', 'numgoodframes');


% load([path 'MyTimestamps.mat'])
% temp = load([path 'MyTimestamps.mat']);
% fn = fieldnames(temp);
% antenna_time = temp.(fn{1});
% antenna_time = posixtime(antenna_time)+t_offset;

% M = readmatrix([path 'log.txt']);
% antenna_time = M(:, 2) + t_offset;


load([path 'AF_ITP_results.mat'])

antenna_time = 0:0.955:0.955*(size(AF_ITP_results)-1);


fullpath = fullfile(path, csv);
T = readtable(fullpath);
data = struct();
for col = 1:width(T)
    data.(T.Properties.VariableNames{col}) = T{:, col};
end

data.datetime_utc_.TimeZone = 'UTC';
t1 = data.datetime_utc_(1);
idxsame = find(data.datetime_utc_(1:10) == t1);
offset = (10-length(idxsame));
data.time_millisecond_ = data.time_millisecond_ - data.time_millisecond_(1);
data.UTC_seconds = data.datetime_utc_(1) + seconds(data.time_millisecond_ + offset*100)/1000;
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


angle_offset = mod(90 - ref_direction, 360); %rad2deg(atan2(data.y(1), data.x(1)));
AF_ITP_results(:, 1) = mod(AF_ITP_results(:, 1) + angle_offset +180, 360) - 180;
dt = data.UTC_seconds - data.UTC_seconds(1);

startk = 220;
endk = 680%length(data.x);

antenna_time = antenna_time + data.UTC_seconds(220) + t_offset;

%% Setup video writer
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

figure;
set(gcf, 'Position', [100, 100, 1400, 700]);

% Pre-compute limits
margin = 10;
xlimVals = [min([min(data.x), min(data.y)])-margin, max([max(data.x), max(data.y)])+margin];
ylimVals = [min([min(data.x), min(data.y)])-margin, max([max(data.x), max(data.y)])+margin];
zlimVals = [min([min(data.x), min(data.y)])-margin, max([max(data.x), max(data.y)])+margin];

% Create subplots
subplot(1,2,1);
axXY = gca;
xlabel('X (m)'); ylabel('Y (m)');
title('XY View');
grid on; axis equal;
xlim(xlimVals); ylim(ylimVals);
hold on;
surf(axXY, sX, sY, sZ+antenna_height);

subplot(1,2,2);
axXZ = gca;
xlabel('X (m)'); ylabel('Z (m)');
title('XZ View');
grid on; axis equal;
xlim(xlimVals); ylim(zlimVals);
hold on;
surf(axXZ, sX, sY+antenna_height, sZ);  % Shift Y to Z visually


h_markerXY = plot(axXY, nan, nan, 'xk', 'LineWidth', 2, 'MarkerSize', 10);
h_markerXZ = plot(axXZ, nan, nan, 'xk', 'LineWidth', 2, 'MarkerSize', 10);

h_beamXY = plot(axXY, nan, nan, 'r--', 'LineWidth', 2);
h_beamXZ = plot(axXZ, nan, nan, 'r--', 'LineWidth', 2);

% Persistent paths (GPS vs antenna)
h_pathXY_gps = plot(axXY, nan, nan, 'b-', 'LineWidth', 1.5, 'Color', [0 0 1 0.4]); % GPS path
h_pathXZ_gps = plot(axXZ, nan, nan, 'b-', 'LineWidth', 1.5, 'Color', [0 0 1 0.4]);

h_pathXY_ant = plot(axXY, nan, nan, 'r-', 'LineWidth', 1.5, 'Color', [1 0 0 0.7]); % Antenna path
h_pathXZ_ant = plot(axXZ, nan, nan, 'r-', 'LineWidth', 1.5, 'Color', [1 0 0 0.7]);
last_I = NaN;


h_text = text(axXY, nan, nan, '', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);

for k = startk:endk
    % Find matching antenna frame
    diff = antenna_time - data.UTC_seconds(k);
    nonnegative_idx = find(diff>=0);
    [a, Irel] = min(diff(nonnegative_idx));
    I = nonnegative_idx(Irel);
    % [a, I] = min(abs(antenna_time - data.UTC_seconds(k)));
    df_file = I;
    if a<.5
        curr_az = deg2rad(AF_ITP_results(I, 1));
        curr_el = deg2rad(AF_ITP_results(I, 2));
        r = sqrt(data.x(k)^2 + data.y(k)^2 + data.z(k)^2)*1;
        x = r*cos(curr_el)*cos(curr_az);
        y = r*cos(curr_el)*sin(curr_az);
        z = r*sin(curr_el) + antenna_height;
        % Update beam
        set(h_beamXY, 'XData', [0, x], 'YData', [0, y]);
        set(h_beamXZ, 'XData', [0, x], 'YData', [antenna_height, z]);
    else
        curr_az = NaN;
        curr_el = NaN;
        r = sqrt(data.x(k)^2 + data.y(k)^2 + data.z(k)^2)*1;
        x = r*cos(curr_el)*cos(curr_az);
        y = r*cos(curr_el)*sin(curr_az);
        z = r*sin(curr_el) + antenna_height;
    end



    % Update GPS path
    gps_x = data.x(1:k);
    gps_y = data.y(1:k);
    gps_z = data.z(1:k);
    
    set(h_pathXY_gps, 'XData', gps_x, 'YData', gps_y);
    set(h_pathXZ_gps, 'XData', gps_x, 'YData', gps_z);
    
    % Update antenna path only when index changes
    if ~isnan(curr_az) && ~isnan(curr_el) && (I ~= last_I)
        ant_x = get(h_pathXY_ant, 'XData');
        ant_y = get(h_pathXY_ant, 'YData');
        ant_z = get(h_pathXZ_ant, 'YData');
    
        set(h_pathXY_ant, 'XData', [ant_x, x], 'YData', [ant_y, y]);
        set(h_pathXZ_ant, 'XData', [ant_x, x], 'YData', [ant_z, z]);
    
        last_I = I; % update tracker
    end





    % Update drone marker
    set(h_markerXY, 'XData', data.x(k), 'YData', data.y(k));
    set(h_markerXZ, 'XData', data.x(k), 'YData', data.z(k));

    % Altitude label in XY plot
    set(h_text, 'Position', [data.x(k), data.y(k)], ...
        'String', sprintf('Alt: %.1f m', data.z(k)));
    % sgtitle( num2str(k))

    drawnow;
    pause(0.01);
    sgtitle(['Time = ' num2str(dt(k) - dt(startk)) ' s;    DF File: ' num2str(df_file) ';    Drone GPS Tick = ' num2str(k) ';'])

    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
fprintf('Video saved to: %s\n', videoName);




%% Plot Figures
drone_az = rad2deg(atan2(data.y, data.x));
drone_el = rad2deg(atan2(data.z-antenna_height, sqrt(data.x.^2 + data.y.^2)));
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



% figure(12)
% sgtitle(['6 Antenna SNR (Noise Level = ' num2str(noise_level_test) ' dB)'])
% for i = 1:6
%     subplot(2, 3, i)
%     scatter(drone_loc(:, 3), Test_Mag(i, :)-noise_level_test)
% 
%     [~, sorted_idx] = sort(drone_loc(:, 3));
%     drone_distance_sort = drone_loc(sorted_idx, 3);
%     SNR = (Test_Mag(i, :)-noise_level_test)';
%     SNR_sorted = SNR(sorted_idx);
% 
%     hold on
%     X = 1 ./ drone_distance_sort;        % compute 1/r
%     Y = SNR_sorted;
%     validIdx = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
%     X = X(validIdx);
%     Y = Y(validIdx);
% 
% 
%     p = polyfit(X, Y, 1);
%     Y_fit = polyval(p, X);
% 
%     % Compute R^2
%     SS_res = sum((Y - Y_fit).^2);         % Residual sum of squares
%     SS_tot = sum((Y - mean(Y)).^2);       % Total sum of squares
%     R_squared = 1 - (SS_res / SS_tot);
% 
%     plot(1./X, Y_fit)
%     eqnStr = sprintf('Fit: P = %.2f / r + %.3f', p(1), p(2));
%     rsquare = sprintf(['R^{2} = ' num2str(R_squared)]);
%     plot(NaN, NaN, 'w')
%     plot(NaN, NaN, 'w')
% 
%     ylabel('Power')
%     xlabel('Distance (m)')
%     ylim([0 2.5*10^4])
%     title(['Antenna ' num2str(i)])
%     grid on
%     legend('Raw Data', '1/r Fit', eqnStr, rsquare, 'Location', 'best')
% 
% 
%     ylabel('SNR (dB)')
%     xlabel('Distance (m)')
%     ylim([0 55])
%     title(['Antenna ' num2str(i)])
%     grid on
% end
% set(gcf, 'Position', [100, 100, 1400, 700]);

% figure(13)
% sgtitle('6 Antenna Power Level')
% for i = 1:6
%     subplot(2, 3, i)
%     scatter(drone_loc(:, 3), 10.^(Test_Mag(i, :)/20))
% 
%     [~, sorted_idx] = sort(drone_loc(:, 3));
%     drone_distance_sort = drone_loc(sorted_idx, 3);
%     mag_sorted = 10.^(Test_Mag(i, sorted_idx)/20)';
% 
%     hold on
%     X = 1 ./ drone_distance_sort.^2;        % compute 1/r^2
%     Y = mag_sorted;
%     validIdx = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
%     X = X(validIdx);
%     Y = Y(validIdx);
% 
% 
%     p = polyfit(X, Y, 1);
%     Y_fit = polyval(p, X);
% 
%     % Compute R^2
%     SS_res = sum((Y - Y_fit).^2);         % Residual sum of squares
%     SS_tot = sum((Y - mean(Y)).^2);       % Total sum of squares
%     R_squared = 1 - (SS_res / SS_tot);
% 
%     plot(sqrt(1./X), Y_fit)
%     eqnStr = sprintf('Fit: P = %.2f / r^{2} + %.3f', p(1), p(2));
%     rsquare = sprintf(['R^{2} = ' num2str(R_squared)]);
%     plot(NaN, NaN, 'w')
%     plot(NaN, NaN, 'w')
%     ylabel('Power')
%     xlabel('Distance (m)')
%     ylim([0 2.5*10^4])
%     title(['Antenna ' num2str(i)])
%     grid on
%     legend('Raw Data', '1/r^{2} Fit', eqnStr, rsquare, 'Location', 'best')
% end
% set(gcf, 'Position', [100, 100, 1400, 700]);

error(:, 1) = mod(error(:, 1) + 180, 360) - 180;
error(:, 2) = mod(error(:, 2) + 180, 360) - 180;


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
    % saveas(figure(12), fullfile(path, 'SNR.jpeg'));
    % saveas(figure(13), fullfile(path, 'Power.jpeg'));
    saveas(figure(14), fullfile(path, 'Error_vs_Distance.jpeg'));
    saveas(figure(15), fullfile(path, 'Error_vs_GT.jpeg'));
    saveas(figure(16), fullfile(path, 'Path.jpeg'));


end