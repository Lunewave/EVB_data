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


libpath = 'U:\Direction_Finding\20251107_Marana_915MHz_CalibrationLibrary_360AZ_66_-6_EL_step3\data';
%%%%%%%%%%% FIXED PARAMETERS %%%%%%%%%%%%%
AZ_data = AZ_start:AZ_step:AZ_end;
AZ_steps = length(AZ_data);
EL_data = EL_start:EL_step:EL_end;
EL_steps = length(EL_data);
numpeaks2check = 1; %# of peaks to check in each dimension of angle interpolation
%%%%%%%%%%%% LIBRARY %%%%%%%%%%%%%%%%%%%%%
offset = 1;
frequency = 915; %MHz
lib_cache = fullfile(libpath, [num2str(frequency/1000) 'GHz_cached_library_data.mat']);
if isfile(lib_cache)
    load(lib_cache, 'Lib_Mag', 'Lib_Phase', 'Lib_Complex', 'Lib_Metadata');
else
    [Lib_Mag, Lib_Phase, Lib_Complex] = Load_FALCON_EVB_Data_900MHz(libpath, AZ_steps, EL_steps, offset, frequency/1000);
    Lib_Metadata = struct(); Lib_Metadata.Azimuth = [AZ_start AZ_step AZ_end]; Lib_Metadata.Elevation = [EL_start EL_step EL_end]; 
    Lib_Metadata.Frequency = '2.456 GHz'; Lib_Metadata.Dimensions = [6 AZ_steps EL_steps]; Lib_Metadata.DimensionString = {'Antennas', 'Azimuth', 'Elevation'};
    save(lib_cache, 'Lib_Mag', 'Lib_Phase', 'Lib_Complex', 'Lib_Metadata');
end
%%%%%%%%%%%%VNA LIB%%%%%%%%%%%%%%%%%
Datadir = 'U:\Antenna test\05142025_AZ_EL_MaranaTest';
filename = '\S_parameters_SingleSweep_FieldTest_WithLens_AZ_EL_';

S21_mag_data_all = zeros(6, 53, 17, 13);
S21_phase_data_all = zeros(6, 53, 17, 13);

for ant = 1:6
    load([Datadir filename 'antenna' num2str(ant) '.mat']);
    output = reshape(output, 5, 53, 17, 13);
    S21_mag_data_all(ant,:,:,:) = reshape(output(2,:,:,:), [1, 53, 17, 13]);
    S21_phase_data_all(ant,:,:,:) = reshape(output(3,:,:,:), [1, 53, 17, 13]);
end

S21_complex=db2mag(S21_mag_data_all).*exp(j.*S21_phase_data_all/180*pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dir = '.\915MHz\';

csvFiles = dir(fullfile(Dir, '*.csv'));
csvData = struct();
for k = 1:length(csvFiles)
    fileName = erase(csvFiles(k).name, '.csv');  % Remove .csv extension for valid field names
    filePath = fullfile(Dir, csvFiles(k).name);
    csvData.(matlab.lang.makeValidName(fileName)) = readtable(filePath);
end

for i = 1:6
    antenna = sprintf(['A' num2str(i)]);
    E_Complex = csvData.(antenna).re_rETheta__mV_ + 1i*csvData.(antenna).im_rETheta__mV_;
    csvData.(antenna).magnitude = 20*log10(abs(E_Complex));
    csvData.(antenna).phase = angle(E_Complex) * 180/pi;
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dir = '.\HFSS_SimDataAMP_AZ_Lens\';
% 
% csvFiles = dir(fullfile(Dir, '*.csv'));
% csvData_Amp_AZ = struct();
% for k = 1:length(csvFiles)
%     fileName = erase(csvFiles(k).name, '.csv');  % Remove .csv extension for valid field names
%     filePath = fullfile(Dir, csvFiles(k).name);
%     csvData_Amp_AZ.(matlab.lang.makeValidName(fileName)) = readtable(filePath);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dir = '.\HFSS_SimDataPhase_AZ_Lens\';
% 
% csvFiles = dir(fullfile(Dir, '*.csv'));
% csvData_Phase_AZ = struct();
% for k = 1:length(csvFiles)
%     fileName = erase(csvFiles(k).name, '.csv');  % Remove .csv extension for valid field names
%     filePath = fullfile(Dir, csvFiles(k).name);
%     csvData_Phase_AZ.(matlab.lang.makeValidName(fileName)) = readtable(filePath);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dir = '.\HFSS_SimDataAMP_EL_Lens\';
% 
% csvFiles = dir(fullfile(Dir, '*.csv'));
% csvData_Amp_EL = struct();
% for k = 1:length(csvFiles)
%     fileName = erase(csvFiles(k).name, '.csv');  % Remove .csv extension for valid field names
%     filePath = fullfile(Dir, csvFiles(k).name);
%     csvData_Amp_EL.(matlab.lang.makeValidName(fileName)) = readtable(filePath);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dir = '.\HFSS_SimDataPhase_EL_Lens\';
% 
% csvFiles = dir(fullfile(Dir, '*.csv'));
% csvData_Phase_EL = struct();
% for k = 1:length(csvFiles)
%     fileName = erase(csvFiles(k).name, '.csv');  % Remove .csv extension for valid field names
%     filePath = fullfile(Dir, csvFiles(k).name);
%     csvData_Phase_EL.(matlab.lang.makeValidName(fileName)) = readtable(filePath);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% freq_indxs_AZ = 181:-1:1;
% freq_indxs_EL = 46:-1:13;

freq_indxs_AZ = 8146:8326;
freq_indxs_EL = 8236:-181:2263;

%% Plot Antenna Patterns
for antenna_ind=1:6
    if AZ_steps>1  

        fieldName = ['A' num2str(antenna_ind)];
        curAntennaAmp = csvData.(fieldName).magnitude;
        curAntennaPhase = csvData.(fieldName).phase;

        LibMag = squeeze(Lib_Mag(antenna_ind,:, 23));
        VNAMag = squeeze(S21_mag_data_all(antenna_ind,2,:, end));

        % Convert to linear power scale
        LibPowerLinear  = 10.^((LibMag - noise_level_cal)/10);
        
        % Compute mean in linear scale
        aSNRLib  = 10 * log10(mean(LibPowerLinear));

        
        figure(100+antenna_ind)
        subplot(1,2,1)
        plot(AZ_data, LibMag);
        grid on; hold on;
        plot([AZ_start:-2:AZ_end], flip(curAntennaAmp(freq_indxs_AZ)));
        % plot([-40:5:40], VNAMag)
        ylim([40 110])
        xlabel('Azimuth Angle (deg)'); ylabel('Magnitude (dB)');
        title(['Antenna ' int2str(antenna_ind) ' magnitude EL = 0']);
        legend(lib_location, 'Simulation Gain Result 915 MHz', 'VNA Data 900 MHz', 'Location', 'best');
        snr_text = sprintf(['Average ' lib_location ' Azimuth SNR: %.2f dB     '], aSNRLib);
        
        uicontrol('Style', 'text', ...
                  'String', snr_text, ...
                  'Units', 'normalized', ...
                  'Position', [0, 0, 1, 0.03], ...  % bottom strip
                  'HorizontalAlignment', 'center', ...
                  'FontSize', 10, ...
                  'BackgroundColor', get(gcf, 'Color'), ...
                  'ForegroundColor', 'k', ...
                  'Tag', 'snr_footer');

        simphase =  unwrap(squeeze(curAntennaPhase(freq_indxs_AZ)), 180) - unwrap(squeeze(csvData.A1.phase(freq_indxs_AZ)), 180);
        dataphase = unwrap(squeeze(Lib_Phase(antenna_ind, :, 23)),180) +shifts1(antenna_ind);
        subplot(1,2,2)
        plot(AZ_data, dataphase - dataphase(end));
        grid on;hold on;
        plot([AZ_start:-2:AZ_end], simphase - simphase(end));
        % plot([-40:5:40], unwrap(squeeze(S21_phase_data_all(antenna_ind, 2, :, end)-S21_phase_data_all(1, 2, :, end)), 180))
        ylim([-240 240])
        xlabel('Angle (deg)');ylabel('Phase (deg)');
        title (['Antenna ' int2str(antenna_ind) ' phase - Antenna 1 phase EL = 0']);
        legend(lib_location, 'Simulation Gain Result 915 MHz', 'VNA Data 900 MHz', 'Location', 'best');
        set(gcf, 'Position', [100, 100, 1200, 500]);
    end

    if EL_steps>1

        fieldName = ['A' num2str(antenna_ind)];
        curAntennaAmp = csvData.(fieldName).magnitude;
        curAntennaPhase = csvData.(fieldName).phase;

        LibMag = squeeze(Lib_Mag(antenna_ind,(length(AZ_data)+1)/2, :));
        VNAMag = squeeze(S21_mag_data_all(antenna_ind,2,7, :));

        % Convert to linear power scale
        LibPowerLinear  = 10.^((LibMag - noise_level_cal)/10);
        
        % Compute mean in linear scale
        aSNRLib  = 10 * log10(mean(LibPowerLinear));

        
        figure(110+antenna_ind)
        subplot(1,2,1)
        plot(EL_data, LibMag);
        grid on; hold on;
        plot([0:2:EL_start], curAntennaAmp(freq_indxs_EL));
        % plot([60:-5:0], VNAMag)
        ylim([40 110])
        xlabel('Elevation Angle (deg)'); ylabel('Magnitude (dB)');
        title(['Antenna ' int2str(antenna_ind) ' magnitude AZ = 0']);
        legend(lib_location, 'Simulation Gain Result 915 MHz', 'VNA Data 900 MHz', 'Location', 'best');
        snr_text = sprintf(['Average ' lib_location ' Elevation SNR: %.2f dB     '], aSNRLib);
        
        uicontrol('Style', 'text', ...
                  'String', snr_text, ...
                  'Units', 'normalized', ...
                  'Position', [0, 0, 1, 0.03], ...  % bottom strip
                  'HorizontalAlignment', 'center', ...
                  'FontSize', 10, ...
                  'BackgroundColor', get(gcf, 'Color'), ...
                  'ForegroundColor', 'k', ...
                  'Tag', 'snr_footer');
        simphase = unwrap(squeeze(curAntennaPhase(freq_indxs_EL)), 180) - unwrap(squeeze(csvData.A1.phase(freq_indxs_EL)), 180);
        dataphase = unwrap(squeeze(Lib_Phase(antenna_ind,(length(AZ_data)+1)/2, :)),180)+shifts11(antenna_ind);
        subplot(1,2,2)
        plot(EL_data, dataphase - dataphase(23));
        grid on;hold on;
        plot([0:2:EL_start], simphase - simphase(1));
        % plot([60:-5:0], unwrap(squeeze(S21_phase_data_all(antenna_ind, 2, 7, :)-S21_phase_data_all(1, 2, 7, :)), 180))
        xlabel('Angle (deg)');ylabel('Phase (deg)');
        ylim([-240 240])
        title (['Antenna ' int2str(antenna_ind) ' phase - Antenna 1 phase AZ = 0']);
        legend(lib_location, 'Simulation Gain Result 915 MHz', 'VNA Data 900 MHz', 'Location', 'best');
        set(gcf, 'Position', [100, 100, 1200, 500]);
    end
end




%% Save Figures
if save_figs
    folderName = [num2str(frequency),'_MHz'];
    newFolderPath = fullfile(libpath, folderName);
    
    if ~exist(newFolderPath, 'dir')
        mkdir(newFolderPath);
    else
        fprintf('Folder "%s" already exists.\n', folderName);
    end
    
    
    if AZ_steps>1
        saveas(figure(101), fullfile(newFolderPath, 'Antenna_1_AZ.jpeg'));
        saveas(figure(102), fullfile(newFolderPath, 'Antenna_2_AZ.jpeg'));
        saveas(figure(103), fullfile(newFolderPath, 'Antenna_3_AZ.jpeg'));
        saveas(figure(104), fullfile(newFolderPath, 'Antenna_4_AZ.jpeg'));
        saveas(figure(105), fullfile(newFolderPath, 'Antenna_5_AZ.jpeg'));
        saveas(figure(106), fullfile(newFolderPath, 'Antenna_6_AZ.jpeg'));
    end
    if EL_steps>1
        saveas(figure(111), fullfile(newFolderPath, 'Antenna_1_EL.jpeg'));
        saveas(figure(112), fullfile(newFolderPath, 'Antenna_2_EL.jpeg'));
        saveas(figure(113), fullfile(newFolderPath, 'Antenna_3_EL.jpeg'));
        saveas(figure(114), fullfile(newFolderPath, 'Antenna_4_EL.jpeg'));
        saveas(figure(115), fullfile(newFolderPath, 'Antenna_5_EL.jpeg'));
        saveas(figure(116), fullfile(newFolderPath, 'Antenna_6_EL.jpeg'));
    end
end


test = NaN(6, 121*25);

for i = 1:6
    data = squeeze(Lib_Mag(i, :, :));
    test(i, :) = data(:);
end


figure()
plot(test')
xlim([0 25*121])
xlabel('File')
ylabel('Magnitude (dB)')