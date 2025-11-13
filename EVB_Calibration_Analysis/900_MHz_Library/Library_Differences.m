close all; clear all; clc;

AZ_start = 180; AZ_end = -180; AZ_step = -3;
EL_start = 66; EL_end = -6; EL_step = -3;
save_figs = 1;
lib_location = 'Calibration Library 915 MHz';
test_location = 'Test Calibration Library 915 MHz';

noise_level_cal = 50;

shifts1 = [0 0 0 0 0 0];
shifts2 = [0 0 0 0 0 0];

shifts11 = [0 0 0 0 0 0];
shifts12 = [0 0 0 0 0 0];


libpath = 'U:\Direction_Finding\20251107_Marana_915MHz_CalibrationLibrary_360AZ_66_-6_EL_step3\data';
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
lib_cache = fullfile(libpath, [num2str(frequency/1000) 'GHz_cached_library_data.mat']);
if isfile(lib_cache)
    load(lib_cache, 'Lib_Mag', 'Lib_Phase', 'Lib_Complex', 'Lib_Metadata');
else
    [Lib_Mag, Lib_Phase, Lib_Complex] = Load_FALCON_EVB_Data_900MHz(libpath, AZ_steps, EL_steps, offset, frequency/1000);
    Lib_Metadata = struct(); Lib_Metadata.Azimuth = [AZ_start AZ_step AZ_end]; Lib_Metadata.Elevation = [EL_start EL_step EL_end]; 
    Lib_Metadata.Frequency = '2.456 GHz'; Lib_Metadata.Dimensions = [6 AZ_steps EL_steps]; Lib_Metadata.DimensionString = {'Antennas', 'Azimuth', 'Elevation'};
    save(lib_cache, 'Lib_Mag', 'Lib_Phase', 'Lib_Complex', 'Lib_Metadata');
end

%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%
offset = 1;
frequency = 915; %MHz
test_cache = fullfile(testpath, [num2str(frequency/1000) 'GHz_cached_test_data.mat']);
if isfile(test_cache)
    load(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex');
else
    [Test_Mag, Test_Phase, Test_Complex] = Load_FALCON_EVB_Data_900MHz(testpath, AZ_steps, EL_steps, offset, frequency/1000);
    save(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex');
end

for i = 1:6
    figure(i)
    t = tiledlayout(1, 2, 'Padding', 'tight', 'TileSpacing', 'tight'); % Tight spacing

    % Compute global min/max for this antenna
    data1 = abs(squeeze(Lib_Phase(i, :, :))).';
    data2 = abs(squeeze(Test_Phase(i, :, :))).';
    cmin = min([data1(:); data2(:)]);
    cmax = max([data1(:); data2(:)]);

    % Calibration Library
    nexttile
    imagesc(AZ_data, EL_data, data1);
    caxis([cmin-5 cmax+5]);
    cb1 = colorbar;
    cb1.Label.String = 'Phase [deg]';
    xlabel('AZ'); ylabel('EL'); 
    title('Calibration Library');

    % Test Library
    nexttile
    imagesc(AZ_data, EL_data, data2);
    caxis([cmin-5 cmax+5]);
    cb2 = colorbar;
    cb2.Label.String = 'Phase [deg]';
    xlabel('AZ'); ylabel('EL'); 
    title('Test Library');

    set(gcf, 'Position', [100, 100, 1400, 700]);
    sgtitle(t, ['Antenna ' num2str(i) ' Phase Comparison'])
end

for i = 1:6
    figure(i+10)
    t = tiledlayout(1, 2, 'Padding', 'tight', 'TileSpacing', 'tight'); % Tight spacing

    % Compute global min/max for this antenna
    data1 = abs(squeeze(Lib_Mag(i, :, :))).';
    data2 = abs(squeeze(Test_Mag(i, :, :))).';
    cmin = min([data1(:); data2(:)]);
    cmax = max([data1(:); data2(:)]);

    % Calibration Library
    nexttile
    imagesc(AZ_data, EL_data, data1);
    caxis([cmin-5 cmax+5]);
    cb1 = colorbar;
    cb1.Label.String = 'Magnitude [dB]';
    xlabel('AZ'); ylabel('EL'); 
    title('Calibration Library');

    % Test Library
    nexttile
    imagesc(AZ_data, EL_data, data2);
    caxis([cmin-5 cmax+5]);
    cb2 = colorbar;
    cb2.Label.String = 'Magnitude [dB]';
    xlabel('AZ'); ylabel('EL'); 
    title('Test Library');

    set(gcf, 'Position', [100, 100, 1400, 700]);
    sgtitle(t, ['Antenna ' num2str(i) ' Magnitude Comparison'])
end


for i = 1:6

    figure(i+20)
    t = tiledlayout(1, 1, 'Padding', 'tight', 'TileSpacing', 'tight'); % Tight spacing

    % Compute global min/max for this antenna
    data1 = abs(squeeze(Lib_Phase(i, :, :))).';
    data2 = abs(squeeze(Test_Phase(i, :, :))).';
    data = abs(data2 - data1);
    cmin = min(data(:));
    cmax = max(data(:));

    % Calibration Library
    nexttile
    imagesc(AZ_data, EL_data, data);
    caxis([cmin cmax+1]);
    cb1 = colorbar;
    cb1.Label.String = 'Phase [deg]';
    xlabel('AZ'); ylabel('EL'); 
    title(['Antenna ' num2str(i) ' Test Library - Calibration Library (abs(Phase))']);

    set(gcf, 'Position', [100, 100, 1400, 700]);


end


for i = 1:6

    figure(i+30)
    t = tiledlayout(1, 1, 'Padding', 'tight', 'TileSpacing', 'tight'); % Tight spacing

    % Compute global min/max for this antenna
    data1 = abs(squeeze(Lib_Mag(i, :, :))).';
    data2 = abs(squeeze(Test_Mag(i, :, :))).';
    data = abs(data2 - data1);
    cmin = min(data(:));
    cmax = max(data(:));

    % Calibration Library
    nexttile
    imagesc(AZ_data, EL_data, data);
    caxis([cmin cmax+1]);
    cb1 = colorbar;
    cb1.Label.String = 'Magnitude [dB]';
    xlabel('AZ'); ylabel('EL'); 
    title(['Antenna ' num2str(i) ' Test Library - Calibration Library (abs(Magnitude))']);

    set(gcf, 'Position', [100, 100, 1400, 700]);


end

for i = 1:6
    figure(i+40)
    t = tiledlayout(1, 2, 'Padding', 'tight', 'TileSpacing', 'tight'); % Tight spacing

    % Compute global min/max for this antenna
    data1 = abs(squeeze(Lib_Mag(i, :, :))).' - noise_level_cal;
    data2 = abs(squeeze(Test_Mag(i, :, :))).' - noise_level_cal;
    cmin = min([data1(:); data2(:)]);
    cmax = max([data1(:); data2(:)]);

    % Calibration Library
    nexttile
    imagesc(AZ_data, EL_data, data1);
    caxis([cmin-5 cmax+5]);
    cb1 = colorbar;
    cb1.Label.String = 'SNR [dB]';
    xlabel('AZ'); ylabel('EL'); 
    title('Calibration Library');

    % Test Library
    nexttile
    imagesc(AZ_data, EL_data, data2);
    caxis([cmin-5 cmax+5]);
    cb2 = colorbar;
    cb2.Label.String = 'SNR [dB]';
    xlabel('AZ'); ylabel('EL'); 
    title('Test Library');

    set(gcf, 'Position', [100, 100, 1400, 700]);
    sgtitle(t, ['Antenna ' num2str(i) ' SNR Comparison'])
end







%% Plot Antenna Patterns
for antenna_ind=1:6
    if AZ_steps>1  


        LibMag = squeeze(Lib_Mag(antenna_ind,:, 9));
        TestMag = squeeze(Test_Mag(antenna_ind,:, 9));


        % Convert to linear power scale
        LibPowerLinear  = 10.^((LibMag - noise_level_cal)/10);
        TestPowerLinear  = 10.^((TestMag - noise_level_cal)/10);
        
        % Compute mean in linear scale
        aSNRLib  = 10 * log10(mean(LibPowerLinear));
        aSNRTest  = 10 * log10(mean(TestPowerLinear));

        
        figure(100+antenna_ind)
        subplot(1,2,1)
        plot(AZ_data, LibMag);
        grid on; hold on;
        % plot([-40:5:40], VNAMag)
        ylim([40 110])
        xlabel('Azimuth Angle (deg)'); ylabel('Magnitude (dB)');
        title(['Antenna ' int2str(antenna_ind) ' magnitude EL = 0']);
        legend('Lib', 'Test', 'Simulation Gain Result 915 MHz', 'VNA Data 900 MHz', 'Location', 'best');
        snr_text = sprintf(['Average ' lib_location ' Azimuth SNR: %.2f dB     ' ...
                            'Average ' test_location ' Azimuth SNR: %.2f dB'], aSNRLib, aSNRTest);
        
        uicontrol('Style', 'text', ...
                  'String', snr_text, ...
                  'Units', 'normalized', ...
                  'Position', [0, 0, 1, 0.03], ...  % bottom strip
                  'HorizontalAlignment', 'center', ...
                  'FontSize', 10, ...
                  'BackgroundColor', get(gcf, 'Color'), ...
                  'ForegroundColor', 'k', ...
                  'Tag', 'snr_footer');

        libphase = unwrap(squeeze(Lib_Phase(antenna_ind, :, 9)),180) +shifts1(antenna_ind);
        testphase = unwrap(squeeze(Test_Phase(antenna_ind, :, 9)),180) +shifts1(antenna_ind);

        subplot(1,2,2)
        plot(AZ_data, libphase - libphase(end));
        grid on;hold on;
        plot(AZ_data, testphase - testphase(end));
        ylim([-240 240])
        xlabel('Angle (deg)');ylabel('Phase (deg)');
        title (['Antenna ' int2str(antenna_ind) ' phase - Antenna 1 phase EL = 0']);
        legend('Lib', 'Test', 'Simulation Gain Result 915 MHz', 'VNA Data 900 MHz', 'Location', 'best');
        set(gcf, 'Position', [100, 100, 1200, 500]);
    end
end
