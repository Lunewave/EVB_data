close all; clear all; clc;
files = dir('Result_0p915_lib_MUSIC\*.mat');
addpath('C:\Users\WorkStation2\Documents\GitHub\EVB_data\EVB_Calibration_Analysis\900_MHz_Library\MUSIC\Result_0p915_lib_MUSIC')
AF_ITP_results = zeros(length(files), 2);
for i = 1:length(files)
    i
    load([num2str(i) '_0p915_result.mat'])
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