close all; clear all; clc;

AZ_start = -180; AZ_end = 180; AZ_step = 3;
EL_start = 66; EL_end = 0; EL_step = -3;
save_figs = 0;
frequency = 2456; %MHz


libpath = 'U:\Falcon_Project\20250625_MaranaTest_AZ360_EL66_Step3_withLens_withEVB_2.456GHz_CalibrationLibrary';
testpath = 'U:\Falcon_Project\20250617_LWOfficeTest_AZ360_EL66_Step3_withLens_withEVB_2.456GHz_LibraryTest';
%%%%%%%%%%% FIXED PARAMETERS %%%%%%%%%%%%%
AZ_data = AZ_start:AZ_step:AZ_end;
AZ_steps = length(AZ_data);
EL_data = EL_start:EL_step:EL_end;
EL_steps = length(EL_data);
numpeaks2check = 1; %# of peaks to check in each dimension of angle interpolation
%%%%%%%%%%%% LIBRARY %%%%%%%%%%%%%%%%%%%%%
offset = 0;
lib_cache = fullfile(libpath, 'cached_library_data.mat');
if isfile(lib_cache)
    load(lib_cache, 'Lib_Mag', 'Lib_Phase', 'Lib_Complex');
else
    [Lib_Mag, Lib_Phase, Lib_Complex] = Load_FALCON_EVB_Data(libpath, AZ_steps, EL_steps, offset);
    save(lib_cache, 'Lib_Mag', 'Lib_Phase', 'Lib_Complex');
end
%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%
offset = 0;
test_cache = fullfile(testpath, 'cached_test_data.mat');
if isfile(test_cache)
    load(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex');
else
    [Test_Mag, Test_Phase, Test_Complex] = Load_FALCON_EVB_Data(testpath, AZ_steps, EL_steps, offset);
    save(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex');
end

folderName = [num2str(frequency),'_MHz'];
newFolderPath = fullfile(testpath, folderName);

if ~exist(newFolderPath, 'dir')
    mkdir(newFolderPath);
else
    fprintf('Folder "%s" already exists.\n', folderName);
end


for i = 1:6
lib_slice  = squeeze(Lib_Mag(i, :, :));
test_slice = squeeze(Test_Mag(i, :, :));
min_val = min([min(lib_slice(:)), min(test_slice(:))]);
max_val = max([max(lib_slice(:)), max(test_slice(:))]);
clim_min = floor(min_val / 5) * 5;
clim_max = ceil(max_val / 5) * 5;
figure(30 + i)
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
sgtitle(['Antenna ' num2str(i) ' Magnitude']);
nexttile
imagesc(AZ_data, EL_data, squeeze(Lib_Mag(i, :, :)));
xlabel('AZ'); ylabel('EL');
title('Library Magnitude');
clim([clim_min clim_max])
nexttile
imagesc(AZ_data, EL_data, squeeze(Test_Mag(i, :, :)));
xlabel('AZ'); ylabel('EL');
title('Test Magnitude');
clim([clim_min clim_max])
cb = colorbar;
cb.Layout.Tile = 'south';
set(gcf, 'Position',  [200, 200, 1200, 500]);
filename = sprintf('Antenna_%d_Magnitude.jpeg', i);
saveas(figure(30 + i), fullfile(newFolderPath, filename));

figure(40 + i)
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
sgtitle(['Antenna ' num2str(i) ' Phase']);
nexttile
imagesc(AZ_data, EL_data, squeeze(Lib_Phase(i, :, :)));
xlabel('AZ'); ylabel('EL');
title('Library Phase');
nexttile
imagesc(AZ_data, EL_data, squeeze(Test_Phase(i, :, :)));
xlabel('AZ'); ylabel('EL');
title('Test Phase');
set(gcf, 'Position',  [200, 200, 1200, 500]);
filename_phase = sprintf('Antenna_%d_Phase.jpeg', i);
saveas(figure(40 + i), fullfile(newFolderPath, filename_phase));

end




