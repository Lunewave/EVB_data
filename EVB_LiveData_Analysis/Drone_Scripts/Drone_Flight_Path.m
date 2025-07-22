close all; clear all; clc;

%% Load Data
[csv, path] = uigetfile('U:\Falcon_Project\*.csv', 'Select CSV Flight Record');
fullpath = fullfile(path, csv);
T = readtable(fullpath);
data = struct();
for col = 1:width(T)
    data.(T.Properties.VariableNames{col}) = T{:, col};
end

% data.UTC_seconds = posixtime(data.datetime_utc_);
data.UTC_seconds = data.time_millisecond_/1000 - data.time_millisecond_(1)/1000;
% data.UTC_seconds = data.UTC_seconds - data.UTC_seconds(1);
% data.UTC_seconds(2:end) = data.UTC_seconds(2:end)-1;

data.z = (data.altitude_above_seaLevel_feet_  - data.altitude_above_seaLevel_feet_(1))* 0.3048;
[data.x, data.y, ~] = geodetic2ecef_custom(data.latitude, data.longitude, data.z);
data.x = data.x - data.x(1);
data.y = data.y - data.y(1);


%% Setup video writer
videoName = fullfile(path, 'Drone_Flight_Path.mp4');
v = VideoWriter(videoName, 'MPEG-4');
v.FrameRate = 10;
open(v);

%% Create figure for animation
numPoints = length(data.UTC_seconds);
trailLength = 5;  % last 5 points in trail
baseSize = 36;    % largest marker size

figure;
set(gcf, 'Position', [100, 100, 1400, 700]);
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;
margin = 0.1;
xlim([min(data.x)-margin, max(data.x)+margin]);
ylim([min(data.y)-margin, max(data.y)+margin]);
zlim([min(data.z)-margin, max(data.z)+margin]);
title('3D position with shrinking trail');

for k = 1:numPoints
    cla; % Clear axes
    
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid on;
    xlim([min(data.x)-margin, max(data.x)+margin]);
    ylim([min(data.y)-margin, max(data.y)+margin]);
    zlim([min(data.z)-margin, max(data.z)+margin]);
    view(3)
    trailIdx = (k - trailLength + 1):k;
    trailIdx = trailIdx(trailIdx > 0);
    
    hold on;
    % Plot oldest points first, newest last (largest size on top)
    for i = length(trailIdx):-1:1
        sizeFactor = baseSize / (length(trailIdx) - i + 1);
        scatter3(data.x(trailIdx(i)), data.y(trailIdx(i)), data.z(trailIdx(i)), ...
            sizeFactor, 'b', 'filled');
    end
    plot3(data.x(k), data.y(k), 0, 'x', 'Color', [0 0 0], 'LineWidth', 2, 'MarkerSize', 10);
    v_line = linspace(0, data.z(k), 50);
    plot3(data.x(k)*ones(1, 50), data.y(k)*ones(1, 50), v_line, ':')
    text(data.x(k), data.y(k), data.z(k), ...
     ['Altitude: ' num2str(data.z(k), '%.2f') ' m'], ...
     'HorizontalAlignment', 'left', ...
     'VerticalAlignment', 'bottom', ...
     'FontSize', 10, ...
     'Color', [0.5 0.5 0.5]);  % optional: gray color


    
    title(sprintf('Time = %.3f', data.UTC_seconds(k)));
    drawnow;
    pause(0.01);
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
fprintf('Video saved to: %s\n', videoName);



