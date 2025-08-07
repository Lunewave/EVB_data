clc; clear all; close all;

% Define the variable name for timestamps (string)
timestampVarName = 'DroneTest_FlyOutAndBack_TS';

i = 1;
stopLoop = false;

% Setup Control from Host PC to EVB
evb= serialport('COM8', 115200);
evb.Timeout = 1;


write(evb, 's', 'char');
write(evb,newline, 'char');

GetConfirmFromEVB(evb, 'WaitForChar');
disp('Frame 0 captured')


% Initialize empty datetime array
timeStamps = datetime.empty;
timeStamps.TimeZone = 'UTC';
timeStamps.Format = 'yyyy-MM-dd HH:mm:ss.SSSSSS';

hFig = figure('Name', 'Press "q" to stop', ...
              'KeyPressFcn', @(src, event) setappdata(src, 'stopFlag', strcmp(event.Key, 'q')));
setappdata(hFig, 'stopFlag', false);

baseTime = datetime('now', 'TimeZone', 'UTC');
tic;

while ishandle(hFig) && ~getappdata(hFig, 'stopFlag')
    elapsed = toc;
    timeStamps(i) = baseTime + seconds(elapsed);
    write(evb, 's', 'char');
    write(evb,newline, 'char');
    
    GetConfirmFromEVB(evb, 'WaitForChar');

    disp(['Frame ' num2str(i) ' Captured at ' char(timeStamps(i))]);

    i = i + 1;
%     pause(0.25);
end

close(hFig);
disp('Loop stopped by user.');

% Save the timestamps variable with the name specified in timestampVarName
% Assign the timeStamps array to a variable with that name in a struct
S.(timestampVarName) = timeStamps;

% Save the struct S to a MAT-file
save('MyTimestamps.mat', '-struct', 'S');

disp(['Timestamps saved to MyTimestamps.mat as variable "' timestampVarName '".']);

