close all; clear all; clc;


%% ROTATOR LOCATION

save_figs = 1;
ref_lat = 32.45177;       % North is positive
ref_lon = -111.21107;     % West is negative
ref_direction = 138;       % 0 is north, 90 is east, 180 is south. This is the direction that the 0 degree azimuth antenna is pointing.
noise_level_test = 45;
antenna_height = 2;
t_offset = -10; %seconds

%% Load Data
[csv, path] = uigetfile('U:\Direction_Finding\*.csv', 'Select CSV Flight Record');

load([path 'AF_ITP_results.mat'])
M = readmatrix([path 'Session_Log.txt']);
antenna_time = NaN(1, M(end, 1));
antenna_time(M(:, 1)) = (M(:, 2) - M(1, 2))/1000 + t_offset;

fullpath = fullfile(path, csv);
T = readtable(fullpath);
data = struct();
for col = 1:width(T)
    data.(T.Properties.VariableNames{col}) = T{:, col};
end



%% Process Drone Data
% data.datetime_utc_.TimeZone = 'UTC';
% t1 = data.datetime_utc_(1);
% idxsame = find(data.datetime_utc_(1:10) == t1);
% offset = (10-length(idxsame));
% data.time_millisecond_ = data.time_millisecond_ - data.time_millisecond_(1);
% data.UTC_seconds = data.datetime_utc_(1) + seconds(data.time_millisecond_ + offset*100)/1000;
% data.UTC_seconds.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
% data.UTC_seconds = posixtime(data.UTC_seconds);



ref_alt_feet = 1450;
ref_alt_m = ref_alt_feet * 0.3048; % Convert to meters

% Use your lat/lon and ref_lat/ref_lon and ref_alt_m
data.z = (data.altitude_above_seaLevel_feet_  - data.altitude_above_seaLevel_feet_(1))* 0.3048;
[data.x, data.y, ~] = geodetic2enu_custom(data.latitude, data.longitude, data.z, ref_lat, ref_lon, ref_alt_m);
data.time = data.time_millisecond_ - data.time_millisecond_(1);


good_idx = find(abs(data.x) < 1000);
data.x = data.x(good_idx);
data.y = data.y(good_idx);
data.z = data.z(good_idx);
data.time = data.time(good_idx) / 1000;
% data.UTC_seconds = data.UTC_seconds(good_idx);




%% Offset Angle Finding Results
angle_offset = mod(90 - ref_direction, 360);
AF_ITP_results(:, 1) = mod(AF_ITP_results(:, 1) + angle_offset +180, 360) - 180;
dt = data.time;

startk = 1;
endk = length(data.x);

%% Filter out detections from the remote this way

drone_az = rad2deg(atan2(data.y, data.x));
drone_el = rad2deg(atan2(data.z-antenna_height, sqrt(data.x.^2 + data.y.^2)));
drone_r = sqrt(data.x.^2 + data.y.^2 + data.z.^2);

%dont use detections where drone is below ground plane, we are not
%calibrated for it
for i =1:length(drone_el)
    if drone_el(i)<0
        drone_az(i) = NaN;
        drone_el(i) = NaN;
        drone_r(i) = NaN;
    end
end


for ifile = 1:length(AF_ITP_results(:, 1))
    [a, I] = min(abs(antenna_time(ifile) - data.time));
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


for i = 1:length(AF_ITP_results)
    if abs(error(i, 1)) > 90
        error(i, :) = [NaN NaN];
        AF_ITP_results(i, :) = [NaN NaN];
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



%% Use this figure to adjust time offset

figure(11)
subplot(1, 2, 1)
plot(data.time, drone_az, 'Color', [0, 0, 0]);
hold on
plot(antenna_time, AF_ITP_results(:, 1), 'o-')
xlabel('Seconds')
ylabel('Azimuth Angle (degrees)')
legend('Drone Ground Truth', 'FALCON DF Results', 'location', 'best')
subplot(1, 2, 2)
plot(data.time, drone_el, 'Color', [0, 0, 0]);
hold on
plot(antenna_time, AF_ITP_results(:, 2), 'o-')
xlabel('Seconds')
ylabel('Elevation Angle (degrees)')
legend('Drone Ground Truth', 'FALCON DF Results', 'location', 'best')
set(gcf, 'Position', [100, 100, 1400, 700]);


%% Setup video writer
videoName = fullfile(path, 'Drone_Flight_Path_XY_XZ.mp4');
v = VideoWriter(videoName, 'MPEG-4');
v.FrameRate = 10;
open(v);

%% XY + XZ Video
numPoints = length(data.time);
trailLength = 5;  % last 5 points in trail
baseSize = 36;    % largest marker size

[sX, sY, sZ] = sphere;

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
    diff = antenna_time - data.time(k);
    nonnegative_idx = find(diff>=0);
    [a, Irel] = min(diff(nonnegative_idx));
    I = nonnegative_idx(Irel);
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
for i = 1:length(AF_ITP_results) - 1
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
for i = 1:length(AF_ITP_results)
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
    folderName = [num2str(2.4*1000),'_MHz'];
    path = fullfile(path, folderName);
    
    if ~exist(path, 'dir')
        mkdir(path);
    else
        fprintf('Folder "%s" already exists.\n', folderName);
    end
    saveas(figure(10), fullfile(path, 'Error.jpeg'));
    saveas(figure(11), fullfile(path, 'Full_Flight.jpeg'));
    saveas(figure(14), fullfile(path, 'Error_vs_Distance.jpeg'));
    saveas(figure(15), fullfile(path, 'Error_vs_GT.jpeg'));
    saveas(figure(16), fullfile(path, 'Path.jpeg'));


end

figure(17)
subplot(1, 2, 1)
plot(AF_ITP_results(:, 1), 'o-')
hold on
plot(drone_loc(:, 1), 'Color', [0, 0, 0])
ylabel('Azimuth')
xlabel('File')

subplot(1, 2, 2)
plot(AF_ITP_results(:, 2), 'o-')
hold on
plot(drone_loc(:, 2), 'Color', [0, 0, 0])
ylabel('Elevation')
xlabel('File')