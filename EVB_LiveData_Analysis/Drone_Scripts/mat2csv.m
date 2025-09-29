close all; clear all; clc;


%% ROTATOR LOCATION

save_figs = 1;
data_freq = 2.447; %Frequency of test data signal in GHz
ref_lat = 32.45130;       % North is positive
ref_lon = -111.21116;     % West is negative
% ref_direction = 92;       % 0 is north, 90 is east, 180 is south. This is the direction that the 0 degree azimuth antenna is pointing.
noise_level_test = 45;

t_offset = 2;

%% Load Data
[csv, path] = uigetfile('U:\Falcon_Project\*.csv', 'Select CSV Flight Record');


test_cache = fullfile(path, [num2str(data_freq) 'GHz_cached_test_data.mat']);
load(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex', 'num_files', 'numgoodframes');


load([path 'MyTimestamps.mat'])
% After loading MyTimestamps.mat
temp = load([path 'MyTimestamps.mat']);
fn = fieldnames(temp);
antenna_time = temp.(fn{1});
antenna_time = posixtime(antenna_time);
% Save antenna_time to CSV
writematrix(antenna_time, fullfile(path, 'antenna_time.csv'));
