close all; clear all; clc;


%% ROTATOR LOCATION

save_figs = 1;
data_freq = 2.447; %Frequency of test data signal in GHz
ref_lat = 32.45130;       % North is positive
ref_lon = -111.21116;     % West is negative
% ref_direction = 92;       % 0 is north, 90 is east, 180 is south. This is the direction that the 0 degree azimuth antenna is pointing.
noise_level_test = 45;



%% Load Data
[csv, path] = uigetfile('U:\Falcon_Project\*.csv', 'Select CSV Flight Record');
%[csv, path] = uigetfile('D:\Lunewave\Falcon_antenna\Real_measurement_data\Drone_data\*.csv', 'Select CSV Flight Record');


test_cache = fullfile(path, [num2str(data_freq) 'GHz_cached_test_data.mat']);
load(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex', 'num_files', 'numgoodframes');


load([path 'MyTimestamps.mat'])
temp = load([path 'MyTimestamps.mat']);
fn = fieldnames(temp);
antenna_time = temp.(fn{1});
antenna_time = posixtime(antenna_time);
load([path 'AF_ITP_results.mat'])
load([path 'AF_ITP_result_minh.mat'])
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
% data.UTC_seconds = polyval(p, 1:length(data.UTC_seconds)) - 3.5;        % Evaluate the fitted line at x
data.UTC_seconds = polyval(p, 1:length(data.UTC_seconds)) - 3.7;        % Evaluate the fitted line at x



% angle_offset = rad2deg(atan2(data.y(1), data.x(1)));
angle_offset = -6;
% angle_offset = mod(90 - ref_direction, 360);
AF_ITP_results(:, 1) = mod(AF_ITP_results(:, 1) + angle_offset +180, 360) - 180;
AF_ITP_result_minh(:, 1) = mod(AF_ITP_result_minh(:, 1) + angle_offset +180, 360) - 180;
AF_ITP_result_minh=AF_ITP_result_minh(2:end,:);
% figure()
% plot(data.UTC_seconds)
% hold on
% plot(antenna_time)
% grid on;

% 
% time_offset=1.795;%hill
% el_offset=2;%hill

time_offset=1.3;%figure8
el_offset=2;%figure8
numPoints = length(data.UTC_seconds);
% Setup video writer
videoName = fullfile(path, 'Drone_Flight_Path_XY_XZ_minh.mp4');
v = VideoWriter(videoName, 'MPEG-4');
v.FrameRate = 10;
open(v);

%% XY + XZ Video

trailLength = 5;  % last 5 points in trail
baseSize = 36;    % largest marker size

[sX, sY, sZ] = sphere;
antenna_height = 2;

figure;
set(gcf, 'Position', [100, 100, 1600, 800]);
% set(gcf, 'Position', [100, 100, 1200, 1000]);

% Pre-compute limits
margin = 10;
% xlimVals = [min([min(data.x), min(data.y)])-margin, max([max(data.x), max(data.y)])+margin];
% ylimVals = [min([min(data.x), min(data.y)])-margin, max([max(data.x), max(data.y)])+margin];
% zlimVals = [min([min(data.x), min(data.y)])-margin, max([max(data.x), max(data.y)])+margin];

xlimVals=[-50 90];
ylimVals=[-40 40];
xlimVals2=[-50 90];
zlimVals=[-40 40];

% xlimVals=[-140 110];
% ylimVals=[-140 110];
% xlimVals2=[-40 40];
% zlimVals=[-40 40];


subplot(1,2,1); axXY_Minh = gca; title('XY View');
xlabel('X (m)'); ylabel('Y (m)'); grid on; axis equal; xlim(xlimVals); ylim(ylimVals); hold on;
surf(axXY_Minh, sX, sY, sZ+antenna_height);

subplot(1,2,2); axXZ_Minh = gca; title('XZ View');
xlabel('X (m)'); ylabel('Z (m)'); grid on; axis equal; xlim(xlimVals2); ylim(zlimVals); hold on;
surf(axXZ_Minh, sX, sY+antenna_height, sZ);

trailLength = 5;
baseSize = 36;



h_beamXY_Minh = plot(axXY_Minh, nan, nan, 'b--', 'LineWidth', 2);
h_beamXZ_Minh = plot(axXZ_Minh, nan, nan, 'b--', 'LineWidth', 2);
h_markerXY_Minh = plot(axXY_Minh, nan, nan, 'xk', 'LineWidth', 2);
h_markerXZ_Minh = plot(axXZ_Minh, nan, nan, 'xk', 'LineWidth', 2);
h_text_Minh = text(axXY_Minh, nan, nan, '', 'FontSize', 10);

legend(axXY_Minh, ...
    [h_beamXY_Minh, h_markerXY_Minh], ...
    {'DF', 'GPS'}, ...
    'Location', 'northeastoutside');

legend(axXZ_Minh, ...
    [h_beamXZ_Minh, h_markerXZ_Minh], ...
    {'DF', 'GPS'}, ...
    'Location', 'northeastoutside');

for k = 1:numPoints
    % Match timestamps
    [a1, I1] = min(abs(antenna_time - data.UTC_seconds(k))); % Ethan
    [a2, I2] = min(abs(antenna_time - data.UTC_seconds(k) - time_offset)); % Minh

    r = norm([data.x(k), data.y(k), data.z(k)]) * 1.1;

    % === Ethan Beam ===
    if a1 < 4
        az1 = deg2rad(AF_ITP_results(I1, 1));
        el1 = deg2rad(AF_ITP_results(I1, 2));
        x1 = r*cos(el1)*cos(az1);
        y1 = r*cos(el1)*sin(az1);
        z1 = r*sin(el1) + antenna_height;
    else
        x1 = nan; y1 = nan; z1 = nan;
    end

    % === Minh Beam ===
    if a2 < 4
        az2 = deg2rad(AF_ITP_result_minh(I2, 1));
        el2 = deg2rad(AF_ITP_result_minh(I2, 2)+el_offset);
        x2 = r*cos(el2)*cos(az2);
        y2 = r*cos(el2)*sin(az2);
        z2 = r*sin(el2) + antenna_height;
    else
        x2 = nan; y2 = nan; z2 = nan;
    end



    % Update Minh plots
    set(h_beamXY_Minh, 'XData', [0 x2], 'YData', [0 y2]);
    set(h_beamXZ_Minh, 'XData', [0 x2], 'YData', [antenna_height z2]);
    set(h_markerXY_Minh, 'XData', data.x(k), 'YData', data.y(k));
    set(h_markerXZ_Minh, 'XData', data.x(k), 'YData', data.z(k));
    set(h_text_Minh, 'Position', [data.x(k) data.y(k)], 'String', sprintf('Alt: %.1f m', data.z(k)));

    drawnow;
    pause(0.01);

    % Save to video
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

for ifile = 1:length(AF_ITP_result_minh(:, 1))
    [a, I] = min(abs(antenna_time(ifile) - time_offset -  data.UTC_seconds));
    if a<1 %Check to see if the drone time is close enough to the antenna time for good error.
        error_minh(ifile, :) = [AF_ITP_result_minh(ifile, 1) - drone_az(I), AF_ITP_result_minh(ifile, 2) - drone_el(I)];
        drone_loc(ifile, :) = [drone_az(I) drone_el(I) drone_r(I)];
        drone_rxy(ifile) = sqrt(data.x(I).^2 + data.y(I).^2);
    else
        error_minh(ifile, :) = [NaN NaN];
        drone_loc(ifile, :) = [NaN NaN NaN];
        drone_rxy(ifile) = NaN;
    end
end


% figure(10)
% 
% % Top-left: Azimuth Error
% subplot(2, 2, 1);
% plot(abs(error(:, 1)), 'o-')
% xlabel('File')
% ylabel('Azimuth Absolute Error (degrees)')
% title('Azimuth Absolute Error vs File (Ethan)')
% ylim([-3 130])
% grid on
% 
% % Top-right: Elevation Error
% subplot(2, 2, 2);
% plot(abs(error(:, 2)), 'o-')
% xlabel('File')
% ylabel('Elevation Absolute Error (degrees)')
% title('Elevation Absolute Error vs File (Ethan)')
% ylim([-3 25])
% grid on
% % Bottom-left: Azimuth Error (minh)
% subplot(2, 2, 3);
% plot(abs(error_minh(:, 1)), 'o-')
% xlabel('File')
% ylabel('Azimuth Absolute Error (degrees)')
% title('Azimuth Absolute Error vs File (Minh)')
% ylim([-3 130])
% grid on
% % Bottom-right: Elevation Error (minh)
% subplot(2, 2, 4);
% plot(abs(error_minh(:, 2)), 'o-')
% xlabel('File')
% ylabel('Elevation Absolute Error (degrees)')
% title('Elevation Absolute Error vs File (Minh)')
% ylim([-3 25])
% grid on
% % Set figure size
% set(gcf, 'Position', [100, 100, 1400, 800]);
% % 
% % 
figure(11)

% % Top-left: Azimuth - Original Results
% subplot(2, 2, 1)
% plot(data.UTC_seconds - data.UTC_seconds(1), drone_az, 'Color', [0, 0, 0]);
% hold on
% plot(antenna_time - data.UTC_seconds(1), AF_ITP_results(:, 1), 'o-')
% xlabel('Seconds')
% ylabel('Azimuth Angle (degrees)')
% title('Azimuth vs Time (Ethan)')
% legend('Drone Ground Truth', 'FALCON DF Results', 'location', 'best')
% grid on
% 
% % Top-right: Elevation - Original Results
% subplot(2, 2, 2)
% plot(data.UTC_seconds - data.UTC_seconds(1), drone_el, 'Color', [0, 0, 0]);
% hold on
% plot(antenna_time - data.UTC_seconds(1), AF_ITP_results(:, 2), 'o-')
% xlabel('Seconds')
% ylabel('Elevation Angle (degrees)')
% title('Elevation vs Time (Ethan)')
% legend('Drone Ground Truth', 'FALCON DF Results', 'location', 'best')
% grid on

% Bottom-left: Azimuth - Minh Results
subplot(1, 2, 1)
plot(data.UTC_seconds - data.UTC_seconds(1), drone_az, 'Color', [0, 0, 0]);
hold on
plot(antenna_time - data.UTC_seconds(1) - time_offset, AF_ITP_result_minh(:, 1), 'o-')
xlabel('Seconds')
ylabel('Azimuth Angle (degrees)')
title('Azimuth vs Time')
legend('Drone Ground Truth', 'FALCON DF Results', 'location', 'best')
grid on

% Bottom-right: Elevation - Minh Results
subplot(1, 2, 2)
plot(data.UTC_seconds - data.UTC_seconds(1), drone_el, 'Color', [0, 0, 0]);
hold on
plot(antenna_time - data.UTC_seconds(1) - time_offset, AF_ITP_result_minh(:, 2) + el_offset, 'o-')
xlabel('Seconds')
ylabel('Elevation Angle (degrees)')
title('Elevation vs Time')
legend('Drone Ground Truth', 'FALCON DF Results', 'location', 'best')
grid on

% Set figure size
set(gcf, 'Position', [100, 100, 1400, 700]);





% figure(12)
% sgtitle(['6 Antenna SNR (Noise Level = ' num2str(noise_level_test) ' dB)'])
% for i = 1:6
%     subplot(2, 3, i)
%     scatter(drone_loc(:, 3), Test_Mag(i, :)-noise_level_test)
%     ylabel('SNR (dB)')
%     xlabel('Distance (m)')
%     ylim([0 55])
%     title(['Antenna ' num2str(i)])
%     grid on
% end
% set(gcf, 'Position', [100, 100, 1400, 700]);
% % 
% figure(13)
% sgtitle('6 Antenna Power Level')
% for i = 1:6
%     subplot(2, 3, i)
%     scatter(drone_loc(:, 3), 10.^(Test_Mag(i, :)/20))
%     ylabel('Power')
%     xlabel('Distance (m)')
%     ylim([0 2.5*10^4])
%     title(['Antenna ' num2str(i)])
%     grid on
% end
% % set(gcf, 'Position', [100, 100, 1400, 700]);
% % 
% figure(14)
% sgtitle('Error vs Distance')
% 
% % Top-left: Azimuth Error (Original)
% subplot(2, 2, 1)
% scatter(drone_loc(:, 3), abs(error(:, 1)))
% xlabel('Distance (m)')
% ylabel('Absolute Azimuth Error (degrees)')
% title('Azimuth Error (Ethan)')
% ylim([0 140])
% grid on
% 
% % Top-right: Elevation Error (Original)
% subplot(2, 2, 2)
% scatter(drone_loc(:, 3), abs(error(:, 2)))
% xlabel('Distance (m)')
% ylabel('Absolute Elevation Error (degrees)')
% title('Elevation Error (Ethan)')
% ylim([0 30])
% grid on
% 
% % Bottom-left: Azimuth Error (Minh)
% subplot(2, 2, 3)
% scatter(drone_loc(:, 3), abs(error_minh(:, 1)))
% xlabel('Distance (m)')
% ylabel('Absolute Azimuth Error (degrees)')
% title('Azimuth Error (Minh)')
% ylim([0 140])
% grid on
% 
% % Bottom-right: Elevation Error (Minh)
% subplot(2, 2, 4)
% scatter(drone_loc(:, 3), abs(error_minh(:, 2)))
% xlabel('Distance (m)')
% ylabel('Absolute Elevation Error (degrees)')
% title('Elevation Error (Minh)')
% ylim([0 30])
% grid on
% 
% % Set figure size
% set(gcf, 'Position', [100, 100, 1400, 800]);
% 
% % 
% figure(15)
% sgtitle('Absolute Error vs Ground Truth')
% 
% % Top-left: Azimuth Error (Original)
% subplot(2, 2, 1)
% scatter(drone_loc(:, 1), abs(error(:, 1)))
% xlabel('Ground Truth Drone Azimuth (degrees)')
% ylabel('Azimuth Absolute Error (DF - GT) (degrees)')
% title('Azimuth Absolute Error (Ethan)')
% xlim([-180 180])
% ylim([-3 150])
% grid on
% 
% % Top-right: Elevation Error (Original)
% subplot(2, 2, 2)
% scatter(drone_loc(:, 2), abs(error(:, 2)))
% xlabel('Ground Truth Drone Elevation (degrees)')
% ylabel('Elevation Absolute Error (DF - GT) (degrees)')
% title('Elevation Absolute Error (Ethan)')
% xlim([0 90])
% ylim([-3 25])
% grid on
% 
% % Bottom-left: Azimuth Error (Minh)
% subplot(2, 2, 3)
% scatter(drone_loc(:, 1), abs(error_minh(:, 1)))
% xlabel('Ground Truth Drone Azimuth (degrees)')
% ylabel('Azimuth Absolute Error (DF - GT) (degrees)')
% title('Azimuth Absolute Error (Minh)')
% xlim([-180 180])
% ylim([-3 150])
% grid on
% 
% % Bottom-right: Elevation Error (Minh Adjusted)
% subplot(2, 2, 4)
% scatter(drone_loc(:, 2), abs(error_minh(:, 2)))
% xlabel('Ground Truth Drone Elevation (degrees)')
% ylabel('Elevation AbsoluteError (DF - GT) (degrees)')
% title('Elevation Absolute Error (Minh)')
% xlim([0 90])
% ylim([-3 25])
% grid on
% 
% % Set figure size
% set(gcf, 'Position', [100, 100, 1400, 800]);
% 
% 
% figure(16)
% scatter(drone_rxy'.*cosd(AF_ITP_result_minh(:, 1)), drone_rxy'.*sind(AF_ITP_result_minh(:, 1)), 'filled')
% hold on
% for i = 1:num_files - 1
%     x1 = drone_rxy(i)'*cosd(AF_ITP_result_minh(i, 1));
%     y1 = drone_rxy(i)'*sind(AF_ITP_result_minh(i, 1));
%     x2 = drone_rxy(i+1)'*cosd(AF_ITP_result_minh(i+1, 1));
%     y2 = drone_rxy(i+1)'*sind(AF_ITP_result_minh(i+1, 1));
%     line([x1 x2], [y1 y2], 'Color', [1 0 0 0.3], 'LineWidth', 1);  % Simulated transparency
%     vec = [x2 - x1, y2 - y1];
%     vec = vec / norm(vec);  % Normalize
%     perp = [-vec(2), vec(1)];  % Perpendicular
%     L = 0.4;  % Length of arrowhead
%     W = 0.2;  % Width of arrowhead
%     base = [x2, y2] - L * vec;
%     arrow_x = [x2, base(1) + W*perp(1), base(1) - W*perp(1)];
%     arrow_y = [y2, base(2) + W*perp(2), base(2) - W*perp(2)];
%     patch(arrow_x, arrow_y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
% end
% for i = 1:num_files
%     text(drone_rxy(i)'*cosd(AF_ITP_result_minh(i, 1)), drone_rxy(i)'*sind(AF_ITP_result_minh(i, 1)), num2str(i), ...
%         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8)
% end
% xlabel('X (m)')
% ylabel('Y (m)')
% title('Signal Path')
% grid on
% set(gcf, 'Position', [100, 100, 1400, 700]);

figure(17); clf;
clear xlim ylim
% drone_rxy = x(~isnan(drone_rxy));


% Convert polar to Cartesian
x_coords = drone_rxy' .* cosd(AF_ITP_result_minh(:, 1));
y_coords = drone_rxy' .* sind(AF_ITP_result_minh(:, 1));

% Scatter all points
scatter(x_coords, y_coords, 'filled');
hold on;

% % Highlight first point (start)
% plot(x_coords(3), y_coords(3), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
% 
% % Highlight last point (end)
% plot(x_coords(end), y_coords(end), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
% Highlight and label start and end points
x_start = drone_rxy(3)' * cosd(AF_ITP_result_minh(3, 1));
y_start = drone_rxy(3)' * sind(AF_ITP_result_minh(3, 1));
x_end = drone_rxy(end)' * cosd(AF_ITP_result_minh(end, 1));
y_end = drone_rxy(end)' * sind(AF_ITP_result_minh(end, 1));

% Big red points
plot(x_start, y_start, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(x_end, y_end, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% Labels
text(x_start, y_start, 'Start', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');
text(x_end, y_end, 'End', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');


% Draw lines and arrowheads
for i = 1:num_files - 1
    x1 = x_coords(i);
    y1 = y_coords(i);
    x2 = x_coords(i+1);
    y2 = y_coords(i+1);

    % Line segment
    line([x1 x2], [y1 y2], 'Color', [1 0 0 0.3], 'LineWidth', 1);

    % Arrow direction
    vec = [x2 - x1, y2 - y1];
    if norm(vec) == 0, continue; end
    vec = vec / norm(vec);
    perp = [-vec(2), vec(1)];

    % Arrow size
    L = 1.5; W = 0.75;
    base = [x2, y2] - L * vec;
    arrow_x = [x2, base(1) + W*perp(1), base(1) - W*perp(1)];
    arrow_y = [y2, base(2) + W*perp(2), base(2) - W*perp(2)];

    % Arrowhead
    patch(arrow_x, arrow_y, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end

% Number each point
for i = 1:num_files
    text(x_coords(i), y_coords(i), num2str(i-2), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);
end

% Labels and layout
xlabel('X (m)');
ylabel('Y (m)');
title('Signal Flight Path with Arrows and Endpoints');
axis equal;
grid on;
set(gcf, 'Position', [100, 100, 1000, 1000]);  % Square layout
xlim([-50 90]);
ylim([-40 40]);


%% Save Figures
if save_figs
    folderName = [num2str(data_freq*1000),'_MHz_Minh'];
    path = fullfile(path, folderName);
    
    if ~exist(path, 'dir')
        mkdir(path);
    else
        fprintf('Folder "%s" already exists.\n', folderName);
    end
%     saveas(figure(10), fullfile(path, 'Error.jpeg'));
    saveas(figure(11), fullfile(path, 'Full_Flight.jpeg'));
%     saveas(figure(12), fullfile(path, 'SNR.jpeg'));
%     saveas(figure(13), fullfile(path, 'Power.jpeg'));
%     saveas(figure(14), fullfile(path, 'Error_vs_Distance.jpeg'));
%     saveas(figure(15), fullfile(path, 'Error_vs_GT.jpeg'));
    saveas(figure(17), fullfile(path, 'Flight_Path.jpeg'));

end