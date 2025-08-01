close all; clear all; clc;


%% ROTATOR LOCATION

save_figs = 1;
data_freq = 2.447; %Frequency of test data signal in GHz
ref_lat = 32.451553;       % North is positive
ref_lon = -111.21108;     % West is negative
ref_direction = 92;       % 0 is north, 90 is east, 180 is south. This is the direction that the 0 degree azimuth antenna is pointing.
noise_level_test = 45;



%% Load Data
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



% angle_offset = rad2deg(atan2(data.y(1), data.x(1)));
angle_offset = mod(90 - ref_direction, 360);
AF_ITP_results(:, 1) = mod(AF_ITP_results(:, 1) + angle_offset +180, 360) - 180;




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

% Initialize graphics
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
    set(h_text, 'Position', [data.x(k), data.y(k)], ...
        'String', sprintf('Alt: %.1f m', data.z(k)));

    drawnow;
    pause(0.01);
    
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
fprintf('Video saved to: %s\n', videoName);




drone_az = rad2deg(atan2(data.y, data.x));
drone_el = rad2deg(atan2(data.z, sqrt(data.x.^2 + data.y.^2)));
drone_r = sqrt(data.x.^2 + data.y.^2 + data.z.^2);


for ifile = 1:length(AF_ITP_results(:, 1))
    [a, I] = min(abs(antenna_time(ifile) - data.UTC_seconds));
    if a<1 %Check to see if the drone time is close enough to the antenna time for good error.
        error(ifile, :) = [AF_ITP_results(ifile, 1) - drone_az(I), AF_ITP_results(ifile, 2) - drone_el(I)];
        drone_loc(ifile, :) = [drone_az(I) drone_el(I) drone_r(I)];

    else
        error(ifile, :) = [NaN NaN];
        drone_loc(ifile, :) = [NaN NaN NaN];

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
    ylabel('Power')
    xlabel('Distance (m)')
    ylim([0 2.5*10^4])
    title(['Antenna ' num2str(i)])
    grid on
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

end

