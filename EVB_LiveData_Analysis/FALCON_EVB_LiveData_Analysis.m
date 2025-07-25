close all; clear all; clc;

save_figs = 1;
frequency = 2456; %MHz
test_location = 'Drone Test';
noise_level_test = 45;
testpath = 'U:\Falcon_Project\20250711_MaranaTest_AZ360_EL0_Step5_withLens_withEVB_2.456GHz_DroneTest_r-10_h-2';

%%%%%%%%%%% FIXED PARAMETERS %%%%%%%%%%%%%
AZ_start = -180; AZ_end = 180; AZ_step = 3;
EL_start = 66; EL_end = 0; EL_step = -3;
AZ_data = AZ_start:AZ_step:AZ_end;
AZ_steps = length(AZ_data);
EL_data = EL_start:EL_step:EL_end;
EL_steps = length(EL_data);
numpeaks2check = 1; %# of peaks to check in each dimension of angle interpolation
%%%%%%%%%%%% LIBRARY %%%%%%%%%%%%%%%%%%%%%
lib_location = 'Calibration Library';
noise_level_cal = 45;
libpath = 'U:\Falcon_Project\20250625_MaranaTest_AZ360_EL66_Step3_withLens_withEVB_2.456GHz_CalibrationLibrary';
offset = 0;
lib_cache = fullfile(libpath, [num2str(frequency/1000) 'GHz_cached_library_data.mat']);
if isfile(lib_cache)
    load(lib_cache, 'Lib_Mag', 'Lib_Phase', 'Lib_Complex');
else
    [Lib_Mag, Lib_Phase, Lib_Complex] = Load_FALCON_EVB_Data(libpath, AZ_steps, EL_steps, offset);
    save(lib_cache, 'Lib_Mag', 'Lib_Phase', 'Lib_Complex');
end
%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%
offset = 0;
data_freq = 2.357; %Frequency of test data signal in GHz
test_cache = fullfile(testpath, [num2str(data_freq) 'GHz_cached_test_data.mat']);
if isfile(test_cache)
    load(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex', 'num_files', 'numgoodframes');
else
    [Test_Mag, Test_Phase, Test_Complex, num_files, numgoodframes] = Load_FALCON_EVB_LiveData(testpath, offset, data_freq);
    save(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex', 'num_files', 'numgoodframes');
end

% %% Plot Antenna Patterns
% for antenna_ind=1:6
%     noise_level_cal = 45;
%     noise_level_test = 45;
% 
%     shifts1 = [0 0 0 0 0 0];
%     shifts2 = [0 0 0 0 0 0];
% 
%     shifts11 = [0 0 0 0 0 0];
%     shifts12 = [0 0 0 0 0 0];
%     if AZ_steps>1        
%         LibMag = squeeze(Lib_Mag(antenna_ind,:, end));
%         TestMag = squeeze(Test_Mag(antenna_ind,:, end));
%         % Convert to linear power scale
%         LibPowerLinear  = 10.^((LibMag - noise_level_cal)/10);
%         TestPowerLinear = 10.^((TestMag - noise_level_test)/10);
% 
%         % Compute mean in linear scale
%         aSNRLib  = 10 * log10(mean(LibPowerLinear));
%         aSNRTest = 10 * log10(mean(TestPowerLinear));
% 
% 
%         figure(100+antenna_ind)
%         subplot(1,2,1)
%         plot(AZ_data, LibMag);
%         grid on; hold on;
%         plot(-180:5:180, TestMag);
%         xlabel('Azimuth Angle (deg)'); ylabel('Magnitude (dB)');
%         title(['Antenna ' int2str(antenna_ind) ' magnitude EL = 0']);
%         legend(lib_location, test_location, 'Location', 'best');
%         snr_text = sprintf(['Average ' lib_location ' Azimuth SNR: %.2f dB     ' ...
%                             'Average ' test_location ' Azimuth SNR: %.2f dB'], aSNRLib, aSNRTest);
% 
%         uicontrol('Style', 'text', ...
%                   'String', snr_text, ...
%                   'Units', 'normalized', ...
%                   'Position', [0, 0, 1, 0.03], ...  % bottom strip
%                   'HorizontalAlignment', 'center', ...
%                   'FontSize', 10, ...
%                   'BackgroundColor', get(gcf, 'Color'), ...
%                   'ForegroundColor', 'k', ...
%                   'Tag', 'snr_footer');
%         subplot(1,2,2)
%         plot(AZ_data,unwrap(squeeze(Lib_Phase(antenna_ind, :, end)),180)+shifts1(antenna_ind));
%         grid on;hold on;
%         plot(-180:5:180,unwrap(squeeze(Test_Phase(antenna_ind, :, end)),180)+shifts2(antenna_ind));
%         xlabel('Angle (deg)');ylabel('Phase (deg)');
%         title (['Antenna ' int2str(antenna_ind) ' phase - Antenna 1 phase EL = 0']);
%         legend(lib_location, test_location, 'Location', 'best');
%         set(gcf, 'Position', [100, 100, 1200, 500]);
%     end
% 
%     if EL_steps>1
% 
%         LibMag = squeeze(Lib_Mag(antenna_ind,(length(AZ_data)+1)/2, :));
%         TestMag = squeeze(Test_Mag(antenna_ind,(length(AZ_data)+1)/2, :));
%         % Convert to linear power scale
%         LibPowerLinear  = 10.^((LibMag - noise_level_cal)/10);
%         TestPowerLinear = 10.^((TestMag - noise_level_test)/10);
% 
%         % Compute mean in linear scale
%         aSNRLib  = 10 * log10(mean(LibPowerLinear));
%         aSNRTest = 10 * log10(mean(TestPowerLinear));
% 
% 
%         figure(110+antenna_ind)
%         subplot(1,2,1)
%         plot(EL_data, LibMag);
%         grid on; hold on;
%         plot(EL_data, TestMag);
%         xlabel('Elevation Angle (deg)'); ylabel('Magnitude (dB)');
%         title(['Antenna ' int2str(antenna_ind) ' magnitude AZ = 0']);
%         legend(lib_location, test_location, 'Location', 'best');
%         snr_text = sprintf(['Average ' lib_location ' Elevation SNR: %.2f dB     ' ...
%                             'Average ' test_location ' Elevation SNR: %.2f dB'], aSNRLib, aSNRTest);
% 
%         uicontrol('Style', 'text', ...
%                   'String', snr_text, ...
%                   'Units', 'normalized', ...
%                   'Position', [0, 0, 1, 0.03], ...  % bottom strip
%                   'HorizontalAlignment', 'center', ...
%                   'FontSize', 10, ...
%                   'BackgroundColor', get(gcf, 'Color'), ...
%                   'ForegroundColor', 'k', ...
%                   'Tag', 'snr_footer');
%         subplot(1,2,2)
%         plot(EL_data,unwrap(squeeze(Lib_Phase(antenna_ind,(length(AZ_data)+1)/2, :)),180)+shifts11(antenna_ind));
%         grid on;hold on;
%         plot(EL_data,unwrap(squeeze(Test_Phase(antenna_ind,(length(AZ_data)+1)/2, :)),180)+shifts12(antenna_ind));
%         xlabel('Angle (deg)');ylabel('Phase (deg)');
%         title (['Antenna ' int2str(antenna_ind) ' phase - Antenna 1 phase AZ = 0']);
%         legend(lib_location, test_location, 'Location', 'best');
%         set(gcf, 'Position', [100, 100, 1200, 500]);
%     end
% end


%% Angle Finding with Interpolation
%%%%%%%%%%%%  DF code  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AZ_table=AZ_data;
EL_table=[EL_start:EL_step:EL_end];
length_el=length(EL_table);
length_az=length(AZ_table);
AF_results = zeros(num_files, 2);
AF_ITP_results = zeros(num_files, 2);
for ifile= 1:num_files   %-140:2:140;                      % object AZ angle input   [-140:140]

        %Test vector
        test=squeeze(Test_Complex(:, ifile)); 

        if isnan(test)
            AF_results(ifile,:) = [NaN NaN];
            AF_ITP_results(ifile,:) = [NaN NaN];
        else
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
                if peak_az > 1 && peak_az < length_az %&& obj_az_input > AZ_start && obj_az_input < AZ_end
                    a = cf_reshape(peak_az - 1, peak_el);
                    b = cf_reshape(peak_az, peak_el);
                    c = cf_reshape(peak_az + 1, peak_el);
                    if (a ~= c)
                        delta_az = 0.5 * (a - c) / (a - 2*b + c);
                    end
                end
                % Define Elevation Neighborhood w/o wrapping
                delta_el = 0; % Assume no interpolation
                if peak_el > 1 && peak_el < length_el %&& obj_el_input < EL_start && obj_el_input > EL_end
                    a1 = cf_reshape(peak_az, peak_el - 1);
                    b1 = cf_reshape(peak_az, peak_el);
                    c1 = cf_reshape(peak_az, peak_el + 1);
                    if (a1 ~= c1)
                        delta_el = 0.5 * (a1 - c1) / (a1 - 2*b1 + c1);
                    end
                end
                % Initial uninterpolated peak positions
                AZ_peak(i) = AZ_table(peak_az);
                EL_peak(i) = EL_table(peak_el);
                % Final interpolated peak positions
                AZ_peak_itp(i) = AZ_peak(i) + delta_az * AZ_step;
                EL_peak_itp(i) = EL_peak(i) + delta_el * EL_step;
    
            end
            % Save results
            AF_results(ifile,:) = [AZ_peak(1) EL_peak(1)];
            AF_ITP_results(ifile,:) = [AZ_peak_itp(1) EL_peak_itp(1)];
        end

end

AF_results=mod(AF_results+180,360)-180;
AF_ITP_results=mod(AF_ITP_results+180,360)-180;



figure(1000)
set(gcf, 'Position',  [200, 200, 800, 600]);
if EL_steps == 1
    plot(AZ_data, mag2db(abs(coe(:,1:end).')))
    xlabel('AZ');ylabel('Correlation Mag')
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
title(['Frame:  ' num2str(ifile)]);
%% Plot Figures
figure(1)
sgtitle(sprintf('Azimuth Results'));
subplot(1, 2, 1)
plot(AF_results(:, 1), 'o-')
title('No Interpolation')
xlabel('Frame')
ylabel('Azimuth Angle')
grid on
legend('Test Data', 'Ground Truth', 'location', 'best')
subplot(1, 2, 2)
plot(AF_ITP_results(:, 1), 'o-')
title('Interpolation')
xlabel('Frame')
ylabel('Azimuth Angle')
grid on
legend('Test Data', 'Ground Truth', 'location', 'best')
set(gcf, 'Position', [100, 100, 1400, 700]);


figure(2)
sgtitle(sprintf('Elevation Results'));
subplot(1, 2, 1)
plot(AF_results(:, 2), 'o-')
title('No Interpolation')
xlabel('Frame')
ylabel('Elevation Angle')
grid on
legend('Test Data', 'Ground Truth', 'location', 'best')
subplot(1, 2, 2)
plot(AF_ITP_results(:, 2), 'o-')
title('Interpolation')
xlabel('Frame')
ylabel('Elevation Angle')
grid on
legend('Test Data', 'Ground Truth', 'location', 'best')
set(gcf, 'Position', [100, 100, 1400, 700]);


figure(4)
plot(numgoodframes)
grid on
xlabel('File #')
ylabel('Frames with Strong Signal')
title('Number of Frames with Signal')
set(gcf, 'Position', [100, 100, 1400, 700]);


figure(5)
scatter(AF_ITP_results(:, 1), AF_ITP_results(:, 2), 'filled')
hold on
for i = 1:num_files
    text(AF_ITP_results(i, 1), AF_ITP_results(i, 2), num2str(i), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8)
end

xlabel('Azimuth Angle (deg)')
ylabel('Elevation Angle (deg)')
title('Signal Path')
grid on
set(gcf, 'Position', [100, 100, 1400, 700]);







%% Save Figures
if save_figs
    folderName = [num2str(data_freq*1000),'_MHz'];
    newFolderPath = fullfile(testpath, folderName);
    
    if ~exist(newFolderPath, 'dir')
        mkdir(newFolderPath);
    else
        fprintf('Folder "%s" already exists.\n', folderName);
    end
    saveas(figure(1), fullfile(newFolderPath, 'Azimuth.jpeg'));
    saveas(figure(2), fullfile(newFolderPath, 'Elevation.jpeg'));
    saveas(figure(5), fullfile(newFolderPath, 'Path.jpeg'));
end









% truth = (-180:AZ_step:180);
% 
% a = AF_ITP_results(:, 1).' - truth;
% 
% idx = find(abs(a) < 10);
% 
% a_good = a(idx);
% truth_good = truth(idx);



