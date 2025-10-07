close all; clear all; clc;

AZ_start = 180; AZ_end = -180; AZ_step = -3;
EL_start = 66; EL_end = -6; EL_step = -3;
save_figs = 1;
lib_location = 'Calibration Library 915 MHz';

noise_level_cal = 45;

shifts1 = [0 0 0 0 0 0];
shifts2 = [0 0 0 0 0 0];

shifts11 = [0 0 0 0 0 0];
shifts12 = [0 0 0 0 0 0];


libpath = 'U:\Direction_Finding\20250924_MaranaCalibrationLibrary_915MHz_360AZ_66_to_-6EL';
testpath = 'U:\Direction_Finding\20250930_MaranaCalibrationLibrary_915MHz_360AZ_66_to_-6EL';

%%%%%%%%%%% FIXED PARAMETERS %%%%%%%%%%%%%
AZ_data = AZ_start:AZ_step:AZ_end;
AZ_steps = length(AZ_data);
EL_data = EL_start:EL_step:EL_end;
EL_steps = length(EL_data);
numpeaks2check = 1; %# of peaks to check in each dimension of angle interpolation
%%%%%%%%%%%% LIBRARY %%%%%%%%%%%%%%%%%%%%%
offset = 1;
frequency = 915; %MHz
lib_cache = fullfile(libpath, [num2str(frequency/1000) 'GHz_cached_library_data_+-5f.mat']);
if isfile(lib_cache)
    load(lib_cache, 'Lib_Mag', 'Lib_Phase', 'Lib_Complex', 'Frequencies');
else
    [Lib_Mag, Lib_Phase, Lib_Complex, Frequencies] = Load_FALCON_EVB_Data_900MHz_MultiFrequency(libpath, AZ_steps, EL_steps, offset, frequency/1000);
    save(lib_cache, 'Lib_Mag', 'Lib_Phase', 'Lib_Complex', 'Frequencies');
end
%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%
offset = 1;
frequency = 915; %MHz
test_cache = fullfile(testpath, [num2str(frequency/1000) 'GHz_cached_test_data_+-5f.mat']);
if isfile(test_cache)
    load(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex', 'Frequencies');
else
    [Test_Mag, Test_Phase, Test_Complex] = Load_FALCON_EVB_Data_900MHz_MultiFrequency(testpath, AZ_steps, EL_steps, offset, frequency/1000);
    save(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex', 'Frequencies');
end