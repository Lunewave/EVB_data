close all; clear all; clc;

AZ_start = 180; AZ_end = -180; AZ_step = -3;
EL_start = 66; EL_end = -6; EL_step = -3;
save_figs = 1;
lib_location = 'Calibration Library 915 MHz';
test_location = 'Test Calibration Library 915 MHz';

noise_level_cal = 45;

shifts1 = [0 0 0 0 0 0];
shifts2 = [0 0 0 0 0 0];

shifts11 = [0 0 0 0 0 0];
shifts12 = [0 0 0 0 0 0];


libpath = 'U:\Direction_Finding\20250924_MaranaCalibrationLibrary_915MHz_360AZ_66_to_-6EL';
testpath = 'U:\Direction_Finding\20250930_MaranaCalibrationLibrary_915MHz_360AZ_66_to_-6EL';


EL_idx = 23;
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


freq_indxs_AZ = 8146:8326;
freq_indxs_EL = 8236:-181:2263;

%% Plot Antenna Patterns
for antenna_ind=1:6
    if AZ_steps>1  

        fieldName = ['A' num2str(antenna_ind)];
        curAntennaAmp = csvData.(fieldName).magnitude;
        curAntennaPhase = csvData.(fieldName).phase;

        LibMag = squeeze(Lib_Mag(antenna_ind,:, EL_idx));
        TestMag = squeeze(Test_Mag(antenna_ind,:, EL_idx));

        VNAMag = squeeze(S21_mag_data_all(antenna_ind,2,:, end));

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
        plot(AZ_data, TestMag)
        plot([AZ_start:-2:AZ_end], flip(curAntennaAmp(freq_indxs_AZ)));
        % plot([-40:5:40], VNAMag)
        ylim([40 110])
        xlabel('Azimuth Angle (deg)'); ylabel('Magnitude (dB)');
        title(['Antenna ' int2str(antenna_ind) ' magnitude EL = 0']);
        legend(lib_location, test_location, 'Simulation Gain Result 915 MHz', 'VNA Data 900 MHz', 'Location', 'best');
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

        simphase =  unwrap(squeeze(curAntennaPhase(freq_indxs_AZ)), 180) - unwrap(squeeze(csvData.A1.phase(freq_indxs_AZ)), 180);
        libphase = unwrap(squeeze(Lib_Phase(antenna_ind, :, EL_idx)),180) +shifts1(antenna_ind);
        testphase = unwrap(squeeze(Test_Phase(antenna_ind, :, EL_idx)),180) +shifts1(antenna_ind);

        subplot(1,2,2)
        plot(AZ_data, libphase - libphase(end));
        grid on;hold on;
        plot(AZ_data, testphase - testphase(end));
        plot([AZ_start:-2:AZ_end], simphase - simphase(end));
        % plot([-40:5:40], unwrap(squeeze(S21_phase_data_all(antenna_ind, 2, :, end)-S21_phase_data_all(1, 2, :, end)), 180))
        ylim([-270 270])
        xlabel('Angle (deg)');ylabel('Phase (deg)');
        title (['Antenna ' int2str(antenna_ind) ' phase - Antenna 1 phase EL = 0']);
        legend(lib_location, test_location, 'Simulation Gain Result 915 MHz', 'VNA Data 900 MHz', 'Location', 'best');
        set(gcf, 'Position', [100, 100, 1200, 500]);
    end

    if EL_steps>1

        fieldName = ['A' num2str(antenna_ind)];
        curAntennaAmp = csvData.(fieldName).magnitude;
        curAntennaPhase = csvData.(fieldName).phase;

        LibMag = squeeze(Lib_Mag(antenna_ind,(length(AZ_data)+1)/2, :));
        TestMag = squeeze(Test_Mag(antenna_ind,(length(AZ_data)+1)/2, :));
        VNAMag = squeeze(S21_mag_data_all(antenna_ind,2,7, :));

        % Convert to linear power scale
        LibPowerLinear  = 10.^((LibMag - noise_level_cal)/10);
        TestPowerLinear  = 10.^((TestMag - noise_level_cal)/10);
        
        % Compute mean in linear scale
        aSNRLib  = 10 * log10(mean(LibPowerLinear));
        aSNRTest  = 10 * log10(mean(TestPowerLinear));

        
        figure(110+antenna_ind)
        subplot(1,2,1)
        plot(EL_data, LibMag);
        grid on; hold on;
        plot(EL_data, TestMag)
        plot([0:2:EL_start], curAntennaAmp(freq_indxs_EL));
        % plot([60:-5:0], VNAMag)
        ylim([40 110])
        xlabel('Elevation Angle (deg)'); ylabel('Magnitude (dB)');
        title(['Antenna ' int2str(antenna_ind) ' magnitude AZ = 0']);
        legend(lib_location, test_location, 'Simulation Gain Result 915 MHz', 'VNA Data 900 MHz', 'Location', 'best');
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


        simphase = unwrap(squeeze(curAntennaPhase(freq_indxs_EL)), 180) - unwrap(squeeze(csvData.A1.phase(freq_indxs_EL)), 180);
        libphase = unwrap(squeeze(Lib_Phase(antenna_ind,(length(AZ_data)+1)/2, :)),180)+shifts11(antenna_ind);
        testphase = unwrap(squeeze(Test_Phase(antenna_ind,(length(AZ_data)+1)/2, :)),180)+shifts11(antenna_ind);
        subplot(1,2,2)
        plot(EL_data, libphase - libphase(EL_idx));
        grid on;hold on;
        plot(EL_data, testphase - testphase(EL_idx));
        plot([0:2:EL_start], simphase - simphase(1));
        % plot([60:-5:0], unwrap(squeeze(S21_phase_data_all(antenna_ind, 2, 7, :)-S21_phase_data_all(1, 2, 7, :)), 180))
        xlabel('Angle (deg)');ylabel('Phase (deg)');
        ylim([-240 240])
        title (['Antenna ' int2str(antenna_ind) ' phase - Antenna 1 phase AZ = 0']);
        legend(lib_location, test_location, 'Simulation Gain Result 915 MHz', 'VNA Data 900 MHz', 'Location', 'best');
        set(gcf, 'Position', [100, 100, 1200, 500]);
    end
end
%% Angle Finding with Interpolation
%%%%%%%%%%%%  DF code  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AZ_table=AZ_data;
EL_table=[EL_start:EL_step:EL_end];
length_el=length(EL_table);
length_az=length(AZ_table);
for obj_az_input=  1:length(AZ_table)   %-140:2:140;                      % object AZ angle input   [-140:140]
    for obj_el_input= 1:length(EL_table)    %-6:2:6;                        % object EL angle input   [-6:6]

        obj_az=AZ_table(obj_az_input)
        obj_el=EL_table(obj_el_input)
        %Test vector
        test=squeeze(Test_Complex(:,obj_az_input,obj_el_input)); 


        weighting=[1 1 1 1 1 1];
        for az=1:AZ_steps
            for el=1:EL_steps
                r=corrcoef(test.*weighting,squeeze(Lib_Complex(:,az,el)).*weighting);
                coe(az,el)=r(1,2);
            end
        end
        cf_reshape=reshape(abs(coe),[length_az length_el]);
        % Find peaks and use max peak
        [pks,locs_y,locs_x]=peaks2(abs(coe),'MinPeakHeight',0.0,'MinPeakDistance',10);
        best_error = inf;
        numpeaks2check = min(numpeaks2check, length(pks));
        peaks_az = locs_y(1:numpeaks2check);
        peaks_el = locs_x(1:numpeaks2check);
        az_angs = AZ_step*(peaks_az-1) + AZ_start;
        el_angs = EL_start + EL_step*(peaks_el-1);
        


        for i = 1:numpeaks2check
            peak_az = locs_y(i);
            peak_el = locs_x(i);
             % Define Azimuth Neighborhood w/o wrapping
            delta_az = 0; % Assume no interpolation
            if peak_az > 1 && peak_az < length_az && obj_az < AZ_start && obj_az > AZ_end
                a = cf_reshape(peak_az - 1, peak_el);
                b = cf_reshape(peak_az, peak_el);
                c = cf_reshape(peak_az + 1, peak_el);
                if (a ~= c)
                    delta_az = 0.5 * (a - c) / (a - 2*b + c);
                end
            end
            % Define Elevation Neighborhood w/o wrapping
            delta_el = 0; % Assume no interpolation
            if peak_el > 1 && peak_el < length_el && obj_el < EL_start && obj_el > EL_end
                a1 = cf_reshape(peak_az, peak_el - 1);
                b1 = cf_reshape(peak_az, peak_el);
                c1 = cf_reshape(peak_az, peak_el + 1);
                if (a1 ~= c1)
                    delta_el = 0.5 * (a1 - c1) / (a1 - 2*b1 + c1);
                end
            end
            % Initial uninterpolated peak positions
            AZ_peak = AZ_table(peak_az);
            EL_peak = EL_table(peak_el);
            % Final interpolated peak positions
            AZ_peak_itp = AZ_peak + delta_az * AZ_step;
            EL_peak_itp = EL_peak + delta_el * EL_step;

            az_error_itp = abs(AZ_peak_itp - obj_az);
            el_error_itp = abs(EL_peak_itp - obj_el);
            total_error = az_error_itp + el_error_itp;
            if total_error < best_error
                best_error = total_error;
                best_AZ_peak = AZ_peak;
                best_EL_peak = EL_peak;
                best_AZ_peak_itp = AZ_peak_itp;
                best_EL_peak_itp = EL_peak_itp;
            end
        end
        % Save errors
        AZ_err(obj_az_input,obj_el_input) = best_AZ_peak - obj_az;
        AZ_err_ITP(obj_az_input,obj_el_input) = best_AZ_peak_itp - obj_az;
        EL_err(obj_az_input,obj_el_input) = best_EL_peak - obj_el; 
        EL_err_ITP(obj_az_input,obj_el_input) = best_EL_peak_itp - obj_el;
    end
end

AZ_err=mod(AZ_err+180,360)-180;
AZ_err_ITP=mod(AZ_err_ITP+180,360)-180;
EL_err=mod(EL_err+180,360)-180;
EL_err_ITP=mod(EL_err_ITP+180,360)-180;


figure(1)
set(gcf, 'Position',  [200, 200, 800, 600]);
if EL_steps == 1
    plot(AZ_data, mag2db(abs(coe(:,1:end).')))
    xlabel('AZ');ylabel('Correlation Mag')
    % ylim([-2 -0.25])
elseif AZ_steps == 1
    plot(EL_data, mag2db(abs(coe(:,1:end).')))
    xlabel('EL');ylabel('Correlation Mag')
    ylim([-2 -0.25])
else
    imagesc(AZ_data,EL_data,mag2db(abs(coe(:,1:end).')))
    caxis([-8 0]);
    colorbar;
    xlabel('AZ');ylabel('EL');
end
title(['Object at:   AZ ' num2str(obj_az) 'deg ;    EL ' num2str(obj_el) ' deg   (Lens1 thicker mounting)']);

% obj_az = -40
% obj_el = 5
%% Plot Angle Errors
%%%%%%%%%%%%%%%%%%%  angle error  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp=AZ_err;


AZ_err_ave=mean(abs(tmp(:)));
AZ_err_std=std(abs(tmp(:)));
AZ_err_max=max(abs(tmp(:)));

tmp=AZ_err_ITP;
AZ_err_ave_ITP=mean(abs(tmp(:)));
AZ_err_std_ITP=std(abs(tmp(:)));
AZ_err_max_ITP=max(abs(tmp(:)));

tmp=EL_err;
EL_err_ave=mean(abs(tmp(:)));
EL_err_std=std(abs(tmp(:)));
EL_err_max=max(abs(tmp(:)));

tmp=EL_err_ITP;
EL_err_ave_ITP=mean(abs(tmp(:)));
EL_err_std_ITP=std(abs(tmp(:)));
EL_err_max_ITP=max(abs(tmp(:)));

% ['AZ error average=' num2str(AZ_err_ave) ' deg;   AZ error std=' num2str(AZ_err_std) ' deg;']
% ['EL error average=' num2str(EL_err_ave) ' deg;   EL error std=' num2str(EL_err_std) ' deg;']
step_error = 10;

figure(10)
subplot(1, 2, 1)
if EL_steps>1 && AZ_steps >1
    imagesc(AZ_data,EL_data,abs(AZ_err.'));
    caxis(abs(AZ_step) * [-0, step_error]);
    cb = colorbar;
    set(cb, 'Ticks', (-step_error * 0) : abs(AZ_step) : (step_error * abs(AZ_step)));
    xlabel('AZ');ylabel('EL');

else
    plot(AZ_data,AZ_err.');
    grid on
    yl = ylim;
    yl_extended = [yl(1)-2, yl(2)+2];
    ylim(yl_extended)
    xlabel('AZ');ylabel('Error (deg)');
end
title(['AZ error at ' num2str(frequency)  ' MHz']);
err_text = sprintf(['Average Azimuth Error: %.3f degrees\n' ...
                    'Average Azimuth Error SD: %.3f degrees\n' ...
                    'Max Azimuth Error: %.3f degrees'], AZ_err_ave, AZ_err_std, AZ_err_max);
annotation('textbox', [0.01, 0.025, 0.3, 0.05], ...
    'String', err_text, ...
    'FitBoxToText', 'on', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 9);
subplot(1, 2, 2)
if EL_steps>1 && AZ_steps >1
    imagesc(AZ_data,EL_data,abs(AZ_err_ITP.'));
    caxis(abs(AZ_step) * [-0, step_error]);
    cb = colorbar;
    set(cb, 'Ticks', (-step_error * 0) : abs(AZ_step) : (step_error * abs(AZ_step)));
    xlabel('AZ');ylabel('EL');

else
    plot(AZ_data,AZ_err_ITP.');
    grid on
    ylim(yl_extended)
    xlabel('AZ');ylabel('Error (deg)');
end
title(['Interpolated AZ error at ' num2str(frequency)  ' MHz']);
err_text = sprintf(['Average Interpolated Azimuth Error: %.3f degrees\n' ...
                    'Average Interpolated Azimuth Error SD: %.3f degrees\n' ...
                    'Max Interpolated Azimuth Error: %.3f degrees'], AZ_err_ave_ITP, AZ_err_std_ITP, AZ_err_max_ITP);
annotation('textbox', [0.47, 0.025, 0.3, 0.05], ...
    'String', err_text, ...
    'FitBoxToText', 'on', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 9);
set(gcf, 'Position', [100, 100, 1400, 700]);



figure(11)
subplot(1, 2, 1)
if EL_steps>1 && AZ_steps >1
    imagesc(AZ_data,EL_data,abs(EL_err.'));
    caxis(abs(EL_step) * [-0, step_error]);
    cb = colorbar;
    set(cb, 'Ticks', (step_error * 0) : abs(EL_step) : (-step_error * EL_step));
    xlabel('AZ');ylabel('EL');
else
    plot(EL_data,EL_err.');
    grid on
    yl = ylim;
    yl_extended = [yl(1)-2, yl(2)+2];
    ylim(yl_extended)
    xlabel('EL');ylabel('Error (deg)');
end
title(['EL error at ' num2str(frequency)  ' MHz']);
err_text = sprintf(['Average Elevation Error: %.3f degrees\n' ...
                    'Average Elevation Error SD: %.3f degrees\n' ...
                    'Max Elevation Error: %.3f degrees'], EL_err_ave, EL_err_std, EL_err_max);
annotation('textbox', [0.01, 0.025, 0.3, 0.05], ...
    'String', err_text, ...
    'FitBoxToText', 'on', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 9);
subplot(1, 2, 2)
if EL_steps>1 && AZ_steps >1
    imagesc(AZ_data,EL_data,abs(EL_err_ITP.'));
    caxis(abs(EL_step) * [-0, step_error]);
    cb = colorbar;
    set(cb, 'Ticks', (step_error * 0) : abs(EL_step) : (-step_error * EL_step));
    xlabel('AZ');ylabel('EL');
else
    plot(EL_data,EL_err_ITP.');
    grid on
    ylim(yl_extended)
    xlabel('EL');ylabel('Error (deg)');
end
title(['Interpolated EL error at ' num2str(frequency)  ' MHz']);
err_text = sprintf(['Average Interpolated Elevation Error: %.3f degrees\n' ...
                    'Average Interpolated Elevation Error SD: %.3f degrees\n' ...
                    'Max Interpolated Elevation Error: %.3f degrees'], EL_err_ave_ITP, EL_err_std_ITP, EL_err_max_ITP);
annotation('textbox', [0.47, 0.025, 0.3, 0.05], ...
    'String', err_text, ...
    'FitBoxToText', 'on', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 9);
set(gcf, 'Position', [100, 100, 1400, 700]);










if EL_steps == 1
    figure(12)
    plot(AZ_data,AZ_err.'+(AZ_data),'o-');
    hold on;
    plot(AZ_data,AZ_data.','o-');
    plot(AZ_data,AZ_err_ITP.'+(AZ_data),'o-');
    xlabel('AZ (deg)');ylabel('AZ angle (deg)');
    grid on; ylim([AZ_end-5 AZ_start+5]);xlim([AZ_end-5 AZ_start+5])
    legend('No Interpolation','Real Angle', 'Interpolated Angle', 'Location', 'best')
    title('Azimuth at EL = 0')

    figure(13)
    plot(AZ_data,AZ_err_ITP.');
    grid on
    xlabel('Azimuth')
    ylabel('Error')
    title('Elevation  = 0 Azimuth Error')
    set(gcf, 'Position',  [200, 200, 1200, 500]);

elseif AZ_steps == 1
    figure(12)
    plot([EL_start:EL_step:EL_end],EL_err((size(AZ_err, 1)+1)/2,(1:abs(EL_step):end))+[EL_start:EL_step:EL_end],'o-');
    hold on;
    plot([EL_start:EL_step:EL_end],[EL_start:EL_step:EL_end].','o-');
    plot([EL_start:EL_step:EL_end],EL_err_ITP((size(AZ_err, 1)+1)/2,(1:abs(EL_step):end))+[EL_start:EL_step:EL_end],'o-');
    xlabel('EL (deg)');ylabel('EL angle (deg)');
    grid on; ylim([EL_end-5 EL_start+5]);xlim([EL_end-5 EL_start+5])
    legend('No Interpolation','Real Angle', 'Interpolated Angle', 'Location', 'best')
    title('Elevation at AZ = 0')

    figure(13)
    plot(EL_data,EL_err_ITP((size(AZ_err, 1)+1)/2,(1:abs(EL_step):end)).');
    grid on
    xlabel('Elevation')
    ylabel('Error')
    title('Azimuth  = 0 Elevation Error')
    set(gcf, 'Position',  [200, 200, 1200, 500]);

else
    figure(12)
    subplot(1,2,1)
    plot(AZ_data,AZ_err(:,EL_idx).'+AZ_data,'o-');
    hold on;
    plot(AZ_data,AZ_data.','o-');
    plot(AZ_data,AZ_err_ITP(:,EL_idx).'+AZ_data,'o-');
    xlabel('AZ (deg)');ylabel('AZ angle (deg)');
    grid on; ylim([AZ_end-5 AZ_start+5]);xlim([AZ_end-5 AZ_start+5])
    legend('No Interpolation','Real Angle', 'Interpolated Angle', 'Location', 'best')
    title('Azimuth at EL = 0')
    subplot(1,2,2)
    plot([EL_start:EL_step:EL_end],EL_err(61,:)+[EL_start:EL_step:EL_end],'o-');
    hold on;
    plot([EL_start:EL_step:EL_end],[EL_start:EL_step:EL_end].','o-');
    plot([EL_start:EL_step:EL_end],EL_err_ITP(61,:)+[EL_start:EL_step:EL_end],'o-');
    xlabel('EL (deg)');ylabel('EL angle (deg)');
    grid on; ylim([EL_end-5 EL_start+5]);xlim([EL_end-5 EL_start+5])
    legend('No Interpolation','Real Angle', 'Interpolated Angle', 'Location', 'best')
    title('Elevation at AZ = 0')
    set(gcf, 'Position',  [200, 200, 1200, 500]);

    figure(13)
    subplot(1, 2, 1)
    plot(AZ_data,AZ_err_ITP(:,EL_idx).');
    grid on
    xlabel('Azimuth')
    ylabel('Error')
    title('Elevation  = 0 Azimuth Error')
    subplot(1, 2, 2)
    plot(EL_data,EL_err_ITP(61,:).');
    grid on
    xlabel('Elevation')
    ylabel('Error')
    title('Azimuth  = 0 Elevation Error')
    set(gcf, 'Position',  [200, 200, 1200, 500]);

end




%% Save Figures
if save_figs
    folderName = [num2str(frequency),'_MHz'];
    newFolderPath = fullfile(testpath, folderName);
    
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
        saveas(figure(10), fullfile(newFolderPath, 'Azimuth_Error.jpeg'))
    end
    if EL_steps>1
        saveas(figure(111), fullfile(newFolderPath, 'Antenna_1_EL.jpeg'));
        saveas(figure(112), fullfile(newFolderPath, 'Antenna_2_EL.jpeg'));
        saveas(figure(113), fullfile(newFolderPath, 'Antenna_3_EL.jpeg'));
        saveas(figure(114), fullfile(newFolderPath, 'Antenna_4_EL.jpeg'));
        saveas(figure(115), fullfile(newFolderPath, 'Antenna_5_EL.jpeg'));
        saveas(figure(116), fullfile(newFolderPath, 'Antenna_6_EL.jpeg'));
        saveas(figure(11), fullfile(newFolderPath, 'Elevation_Error.jpeg'))
    end

    saveas(figure(12), fullfile(newFolderPath, 'AZ0_and_EL0.jpeg'))
end




clear a b c d e1
a = reshape(EL_err_ITP ,[1 AZ_steps*EL_steps]);
b = reshape(AZ_err_ITP ,[1 AZ_steps*EL_steps]);
figure()
plot(a)
grid on
xlabel('File #')
ylabel('EL Error')

c = a>6;
d = find(c == 1);
e1 = mod(d-1, AZ_steps)+1;
e2 = floor((d - 1) / AZ_steps) + 1;   % elevation index (increments every AZ_steps)
az_ang_axis = 180:-3:-180;
el_ang_axis = 66:-3:-6;

az_angs = az_ang_axis(e1);
el_angs = el_ang_axis(e2);
% Histogram with 1-wide bins
figure()
% Bin edges centered on measurement angles
bin_width = 3;  % since your AZ step is 3 degrees
edges = (min(az_ang_axis)-0.5) : 1 : (max(az_ang_axis)+0.5);
histogram(flip(az_angs), edges);
xlabel('AZ Angle');
ylabel('Count');
title('Histogram of AZ angle with EL Err > 6');




figure()
% Bin edges centered on measurement angles
bin_width = 3;  % since your AZ step is 3 degrees
edges = (min(el_ang_axis)-0.5) : 1 : (max(el_ang_axis)+0.5);
histogram(flip(el_angs), edges);
xlabel('EL Angle');
ylabel('Count');
title('Histogram of EL angle with EL Err > 6');




figure()
plot(b)
grid on
xlabel('File #')
ylabel('AZ Error')

c = abs(b)>50;
d = find(c == 1);
e1 = mod(d-1, AZ_steps)+1;
e2 = floor((d - 1) / AZ_steps) + 1;   % elevation index (increments every AZ_steps)
az_ang_axis = 180:-3:-180;
el_ang_axis = 66:-3:-6;

az_angs = az_ang_axis(e1);
el_angs = el_ang_axis(e2);
% Histogram with 1-wide bins
figure()
% Bin edges centered on measurement angles
bin_width = 3;  % since your AZ step is 3 degrees
edges = (min(az_ang_axis)-0.5) : 1 : (max(az_ang_axis)+0.5);
histogram(flip(az_angs), edges);
xlabel('AZ Angle');
ylabel('Count');
title('Histogram of AZ angle with AZ abs(Err) > 50');




figure()
% Bin edges centered on measurement angles
bin_width = 3;  % since your AZ step is 3 degrees
edges = (min(el_ang_axis)-0.5) : 1 : (max(el_ang_axis)+0.5);
histogram(flip(el_angs), edges);
xlabel('EL Angle');
ylabel('Count');
title('Histogram of EL angle with AZ abs(Err) > 50');