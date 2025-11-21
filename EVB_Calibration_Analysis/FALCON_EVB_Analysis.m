close all; clear all; clc;

AZ_start = -180; AZ_end = 180; AZ_step = 3;
EL_start = 66; EL_end = 0; EL_step = -3;
save_figs = 0;
lib_location = 'Calibration Library';
test_location = 'Marana Test';

noise_level_cal = 45;
noise_level_test = 45;

shifts1 = [0 0 0 0 0 0];
shifts2 = [0 0 0 0 0 0];

shifts11 = [0 0 0 0 0 0];
shifts12 = [0 0 0 0 0 0];


libpath = 'U:\Falcon_Project\20250625_MaranaTest_AZ360_EL66_Step3_withLens_withEVB_2.456GHz_CalibrationLibrary';
testpath = 'U:\Falcon_Project\20250626_MaranaTest_AZ360_EL66_Step3_withLens_withEVB_2.456GHz_TestData_skipfirsttwo';
%%%%%%%%%%% FIXED PARAMETERS %%%%%%%%%%%%%
AZ_data = AZ_start:AZ_step:AZ_end;
AZ_steps = length(AZ_data);
EL_data = EL_start:EL_step:EL_end;
EL_steps = length(EL_data);
numpeaks2check = 1; %# of peaks to check in each dimension of angle interpolation
%%%%%%%%%%%% LIBRARY %%%%%%%%%%%%%%%%%%%%%
offset = 0;
frequency = 2456; %MHz
lib_cache = fullfile(libpath, [num2str(frequency/1000) 'GHz_cached_library_data.mat']);
if isfile(lib_cache)
    load(lib_cache, 'Lib_Mag', 'Lib_Phase', 'Lib_Complex', 'Lib_Metadata');
else
    [Lib_Mag, Lib_Phase, Lib_Complex] = Load_FALCON_EVB_Data(libpath, AZ_steps, EL_steps, offset, frequency/1000);
    Lib_Metadata = struct(); Lib_Metadata.Azimuth = [AZ_start AZ_step AZ_end]; Lib_Metadata.Elevation = [EL_start EL_step EL_end]; 
    Lib_Metadata.Frequency = '2.456 GHz'; Lib_Metadata.Dimensions = [6 AZ_steps EL_steps]; Lib_Metadata.DimensionString = {'Antennas', 'Azimuth', 'Elevation'};
    save(lib_cache, 'Lib_Mag', 'Lib_Phase', 'Lib_Complex', 'Lib_Metadata');
end
%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%
offset = 2;
data_freq = 2456; %MHz
test_cache = fullfile(testpath, [num2str(data_freq/1000) 'GHz_cached_test_data.mat']);
if isfile(test_cache)
    load(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex');
else
    [Test_Mag, Test_Phase, Test_Complex] = Load_FALCON_EVB_Data(testpath, AZ_steps, EL_steps, offset, data_freq/1000);
    save(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex');
end



%% Plot Antenna Patterns
for antenna_ind=1:6
    if AZ_steps>1        
        LibMag = squeeze(Lib_Mag(antenna_ind,:, end));
        TestMag = squeeze(Test_Mag(antenna_ind,:, end));
        % Convert to linear power scale
        LibPowerLinear  = 10.^((LibMag - noise_level_cal)/10);
        TestPowerLinear = 10.^((TestMag - noise_level_test)/10);
        
        % Compute mean in linear scale
        aSNRLib  = 10 * log10(mean(LibPowerLinear));
        aSNRTest = 10 * log10(mean(TestPowerLinear));

        
        figure(100+antenna_ind)
        subplot(1,2,1)
        plot(AZ_data, LibMag);
        grid on; hold on;
        plot(AZ_data, TestMag);
        xlabel('Azimuth Angle (deg)'); ylabel('Magnitude (dB)');
        title(['Antenna ' int2str(antenna_ind) ' magnitude EL = 0']);
        legend(lib_location, test_location, 'Location', 'best');
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
        subplot(1,2,2)
        plot(AZ_data,unwrap(squeeze(Lib_Phase(antenna_ind, :, end)),180)+shifts1(antenna_ind));
        grid on;hold on;
        plot(AZ_data,unwrap(squeeze(Test_Phase(antenna_ind, :, end)),180)+shifts2(antenna_ind));
        xlabel('Angle (deg)');ylabel('Phase (deg)');
        title (['Antenna ' int2str(antenna_ind) ' phase - Antenna 1 phase EL = 0']);
        legend(lib_location, test_location, 'Location', 'best');
        set(gcf, 'Position', [100, 100, 1200, 500]);
    end

    if EL_steps>1

        LibMag = squeeze(Lib_Mag(antenna_ind,(length(AZ_data)+1)/2, :));
        TestMag = squeeze(Test_Mag(antenna_ind,(length(AZ_data)+1)/2, :));
        % Convert to linear power scale
        LibPowerLinear  = 10.^((LibMag - noise_level_cal)/10);
        TestPowerLinear = 10.^((TestMag - noise_level_test)/10);
        
        % Compute mean in linear scale
        aSNRLib  = 10 * log10(mean(LibPowerLinear));
        aSNRTest = 10 * log10(mean(TestPowerLinear));

        
        figure(110+antenna_ind)
        subplot(1,2,1)
        plot(EL_data, LibMag);
        grid on; hold on;
        plot(EL_data, TestMag);
        xlabel('Elevation Angle (deg)'); ylabel('Magnitude (dB)');
        title(['Antenna ' int2str(antenna_ind) ' magnitude AZ = 0']);
        legend(lib_location, test_location, 'Location', 'best');
        snr_text = sprintf(['Average ' lib_location ' Elevation SNR: %.2f dB     ' ...
                            'Average ' test_location ' Elevation SNR: %.2f dB'], aSNRLib, aSNRTest);
        
        uicontrol('Style', 'text', ...
                  'String', snr_text, ...
                  'Units', 'normalized', ...
                  'Position', [0, 0, 1, 0.03], ...  % bottom strip
                  'HorizontalAlignment', 'center', ...
                  'FontSize', 10, ...
                  'BackgroundColor', get(gcf, 'Color'), ...
                  'ForegroundColor', 'k', ...
                  'Tag', 'snr_footer');
        subplot(1,2,2)
        plot(EL_data,unwrap(squeeze(Lib_Phase(antenna_ind,(length(AZ_data)+1)/2, :)),180)+shifts11(antenna_ind));
        grid on;hold on;
        plot(EL_data,unwrap(squeeze(Test_Phase(antenna_ind,(length(AZ_data)+1)/2, :)),180)+shifts12(antenna_ind));
        xlabel('Angle (deg)');ylabel('Phase (deg)');
        title (['Antenna ' int2str(antenna_ind) ' phase - Antenna 1 phase AZ = 0']);
        legend(lib_location, test_location, 'Location', 'best');
        set(gcf, 'Position', [100, 100, 1200, 500]);
    end
end


%% Angle Finding with Interpolation
%%%%%%%%%%%%  DF code  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AZ_table=AZ_data;
EL_table=[EL_start:EL_step:EL_end];
length_el=length(EL_table);
length_az=length(AZ_table);
for obj_az_input=  AZ_start:AZ_step:AZ_end   %-140:2:140;                      % object AZ angle input   [-140:140]
    for obj_el_input= EL_start:EL_step:EL_end    %-6:2:6;                        % object EL angle input   [-6:6]

        obj_az=obj_az_input
        obj_el=obj_el_input
        %Test vector
        test=squeeze(Test_Complex(:,(obj_az-AZ_start)/AZ_step + 1,(EL_start-obj_el)/abs(EL_step) + 1)); 


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
            if peak_az > 1 && peak_az < length_az && obj_az_input > AZ_start && obj_az_input < AZ_end
                a = cf_reshape(peak_az - 1, peak_el);
                b = cf_reshape(peak_az, peak_el);
                c = cf_reshape(peak_az + 1, peak_el);
                if (a ~= c)
                    delta_az = 0.5 * (a - c) / (a - 2*b + c);
                end
            end
            % Define Elevation Neighborhood w/o wrapping
            delta_el = 0; % Assume no interpolation
            if peak_el > 1 && peak_el < length_el && obj_el_input < EL_start && obj_el_input > EL_end
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

            az_error_itp = abs(AZ_peak_itp - obj_az_input);
            el_error_itp = abs(EL_peak_itp - obj_el_input);
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
        AZ_err(obj_az_input + AZ_end + 1, abs(EL_start - obj_el_input) + 1) = best_AZ_peak - obj_az_input;
        AZ_err_ITP(obj_az_input + AZ_end + 1, abs(EL_start - obj_el_input) + 1) = best_AZ_peak_itp - obj_az_input;
        EL_err(obj_az_input + AZ_end + 1, abs(EL_start - obj_el_input) + 1) = best_EL_peak - obj_el_input; 
        EL_err_ITP(obj_az_input + AZ_end + 1, abs(EL_start - obj_el_input) + 1) = best_EL_peak_itp - obj_el_input;
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
tmp=AZ_err(1:AZ_step:end,1:abs(EL_step):end);
AZ_err_ave=mean(abs(tmp(:)));
AZ_err_std=std(abs(tmp(:)));
AZ_err_max=max(abs(tmp(:)));

tmp=AZ_err_ITP(1:AZ_step:end,1:abs(EL_step):end);
AZ_err_ave_ITP=mean(abs(tmp(:)));
AZ_err_std_ITP=std(abs(tmp(:)));
AZ_err_max_ITP=max(abs(tmp(:)));

tmp=EL_err(1:AZ_step:end,1:abs(EL_step):end);
EL_err_ave=mean(abs(tmp(:)));
EL_err_std=std(abs(tmp(:)));
EL_err_max=max(abs(tmp(:)));

tmp=EL_err_ITP(1:AZ_step:end,1:abs(EL_step):end);
EL_err_ave_ITP=mean(abs(tmp(:)));
EL_err_std_ITP=std(abs(tmp(:)));
EL_err_max_ITP=max(abs(tmp(:)));

% ['AZ error average=' num2str(AZ_err_ave) ' deg;   AZ error std=' num2str(AZ_err_std) ' deg;']
% ['EL error average=' num2str(EL_err_ave) ' deg;   EL error std=' num2str(EL_err_std) ' deg;']
step_error = 3;

figure(10)
subplot(1, 2, 1)
if EL_steps>1
    imagesc(AZ_data,EL_data,abs(AZ_err(1:AZ_step:end,1:abs(EL_step):end).'));
    caxis(AZ_step * [-0, step_error]);
    cb = colorbar;
    set(cb, 'Ticks', (-step_error * 0) : AZ_step : (step_error * AZ_step));
    xlabel('AZ');ylabel('EL');

else
    plot(AZ_data,AZ_err(1:AZ_step:end,1:abs(EL_step):end).');
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
if EL_steps>1
    imagesc(AZ_data,EL_data,abs(AZ_err_ITP(1:AZ_step:end,1:abs(EL_step):end).'));
    caxis(AZ_step * [-0, step_error]);
    cb = colorbar;
    set(cb, 'Ticks', (-step_error * 0) : AZ_step : (step_error * AZ_step));
    xlabel('AZ');ylabel('EL');

else
    plot(AZ_data,AZ_err_ITP(1:AZ_step:end,1:abs(EL_step):end).');
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
if AZ_steps>1
    imagesc(AZ_data,EL_data,abs(EL_err(1:AZ_step:end,1:abs(EL_step):end).'));
    caxis(-EL_step * [-0, step_error]);
    cb = colorbar;
    set(cb, 'Ticks', (step_error * 0) : -EL_step : (-step_error * EL_step));
    xlabel('AZ');ylabel('EL');
else
    plot(EL_data,EL_err(1:AZ_step:end,1:abs(EL_step):end).');
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
if AZ_steps>1
    imagesc(AZ_data,EL_data,abs(EL_err_ITP(1:AZ_step:end,1:abs(EL_step):end).'));
    caxis(-EL_step * [-0, step_error]);
    cb = colorbar;
    set(cb, 'Ticks', (step_error * 0) : -EL_step : (-step_error * EL_step));
    xlabel('AZ');ylabel('EL');
else
    plot(EL_data,EL_err_ITP(1:AZ_step:end,1:abs(EL_step):end).');
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
    plot(AZ_data,AZ_err(1:AZ_step:end,end).'+AZ_data,'o-');
    hold on;
    plot(AZ_data,AZ_data.','o-');
    plot(AZ_data,AZ_err_ITP(1:AZ_step:end,end).'+AZ_data,'o-');
    xlabel('AZ (deg)');ylabel('AZ angle (deg)');
    grid on; ylim([AZ_start-5 AZ_end+5]);xlim([AZ_start-5 AZ_end+5])
    legend('No Interpolation','Real Angle', 'Interpolated Angle', 'Location', 'best')
    title('Azimuth at EL = 0')

    figure(13)
    plot(AZ_data,AZ_err_ITP(1:AZ_step:end,end).');
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
    plot(AZ_data,AZ_err(1:AZ_step:end,end).'+AZ_data,'o-');
    hold on;
    plot(AZ_data,AZ_data.','o-');
    plot(AZ_data,AZ_err_ITP(1:AZ_step:end,end).'+AZ_data,'o-');
    xlabel('AZ (deg)');ylabel('AZ angle (deg)');
    grid on; ylim([AZ_start-5 AZ_end+5]);xlim([AZ_start-5 AZ_end+5])
    legend('No Interpolation','Real Angle', 'Interpolated Angle', 'Location', 'best')
    title('Azimuth at EL = 0')
    subplot(1,2,2)
    plot([EL_start:EL_step:EL_end],EL_err((size(AZ_err, 1)+1)/2,(1:abs(EL_step):end))+[EL_start:EL_step:EL_end],'o-');
    hold on;
    plot([EL_start:EL_step:EL_end],[EL_start:EL_step:EL_end].','o-');
    plot([EL_start:EL_step:EL_end],EL_err_ITP((size(AZ_err, 1)+1)/2,(1:abs(EL_step):end))+[EL_start:EL_step:EL_end],'o-');
    xlabel('EL (deg)');ylabel('EL angle (deg)');
    grid on; ylim([EL_end-5 EL_start+5]);xlim([EL_end-5 EL_start+5])
    legend('No Interpolation','Real Angle', 'Interpolated Angle', 'Location', 'best')
    title('Elevation at AZ = 0')
    set(gcf, 'Position',  [200, 200, 1200, 500]);

    figure(13)
    subplot(1, 2, 1)
    plot(AZ_data,AZ_err_ITP(1:AZ_step:end,end).');
    grid on
    xlabel('Azimuth')
    ylabel('Error')
    title('Elevation  = 0 Azimuth Error')
    subplot(1, 2, 2)
    plot(EL_data,EL_err_ITP((size(AZ_err, 1)+1)/2,(1:abs(EL_step):end)).');
    grid on
    xlabel('Elevation')
    ylabel('Error')
    title('Azimuth  = 0 Elevation Error')
    set(gcf, 'Position',  [200, 200, 1200, 500]);

end

%% Save Figures
if save_figs
    folderName = [num2str(data_freq),'_MHz'];
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
        saveas(figure(10), fullfile(newFolderPath, 'Azimuth_Error.jpeg'));
    end
    if EL_steps>1
        saveas(figure(111), fullfile(newFolderPath, 'Antenna_1_EL.jpeg'));
        saveas(figure(112), fullfile(newFolderPath, 'Antenna_2_EL.jpeg'));
        saveas(figure(113), fullfile(newFolderPath, 'Antenna_3_EL.jpeg'));
        saveas(figure(114), fullfile(newFolderPath, 'Antenna_4_EL.jpeg'));
        saveas(figure(115), fullfile(newFolderPath, 'Antenna_5_EL.jpeg'));
        saveas(figure(116), fullfile(newFolderPath, 'Antenna_6_EL.jpeg'));
        saveas(figure(11), fullfile(newFolderPath, 'Elevation_Error.jpeg'));
    end
    saveas(figure(12), fullfile(newFolderPath, 'AZ0_and_EL0.jpeg'));
    
    
    % close all
end