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
data_freq = [2.3537 2.3597] ; % Frequency of test data signal in GHz
test_cache = fullfile(testpath, [num2str(data_freq(1)) '-' num2str(data_freq(2)) 'GHz_cached_test_data.mat']);
if isfile(test_cache)
    load(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex', 'num_files', 'freq_range');
else
    [Test_Mag, Test_Phase, Test_Complex, num_files, freq_range] = Load_FALCON_EVB_LiveData_AllData(testpath, offset, data_freq);
    save(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex', 'num_files', 'freq_range');
end

%% Angle Finding with Interpolation
%%%%%%%%%%%%  DF code  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AZ_table=AZ_data;
EL_table=[EL_start:EL_step:EL_end];
length_el=length(EL_table);
length_az=length(AZ_table);
AF_results = NaN(num_files, 1024, length(freq_range), 2);
AF_ITP_results = NaN(num_files, 1024, length(freq_range), 2);
for ifile = 1:num_files   %-140:2:140;                     % object AZ angle input   [-140:140]
    for frame = 1:1024
        disp(ifile)
        disp(frame)
        for frequency = 1:length(freq_range)
    
            %Test vector
            test=squeeze(Test_Complex(:, ifile, frame, frequency)); 
            if any(isnan(test))
                AF_results(ifile,frame,frequency,:) = [NaN NaN];
                AF_ITP_results(ifile,frame,frequency,:) = [NaN NaN];
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
                AF_results(ifile,frame,frequency,:) = [AZ_peak(1) EL_peak(1)];
                AF_ITP_results(ifile,frame,frequency,:) = [AZ_peak_itp(1) EL_peak_itp(1)];
            end
        end
    end
end

AF_results=mod(AF_results+180,360)-180;
AF_ITP_results=mod(AF_ITP_results+180,360)-180;


AF_results_freq_median = squeeze(median(AF_results, 3, 'omitnan'));
AF_ITP_results_freq_median = squeeze(median(AF_ITP_results, 3, 'omitnan'));



figure(1)
hold on
for i=1:73
    scatter(i*ones(1, 1024), squeeze(AF_results_freq_median(i, :, 1)));
end
xlabel('File Number')
ylabel('Azimuth Angle (deg)')
grid on
set(gcf, 'Position', [100, 100, 1400, 700]);





% num_files = 27;
% frames = 1:16:(1024-16); % frames per file
% num_frames = length(frames);
% 
% % Preallocate Y for all files and frames
% Y = zeros(num_files * num_frames, 1);
% 
% for f = 1:num_files
%     % Extract your data from AF_ITP_results_freq_median for each file and frames
%     % Replace the frequency index (3rd dim) as needed, here assumed 1
%     Y(((f-1)*num_frames + 1):(f*num_frames)) = squeeze(AF_ITP_results_freq_median(f, frames, 2));
% end
% 
% % Build X axis with fractional labels per file and frame
% X = [];
% for f = 1:num_files
%     X = [X; f + frames'/1000];  % e.g. 1.001, 1.065, etc.
% end
% 
% % Plot
% figure;
% plot(X, Y, 'o-');
% xlabel('File.Frame Index');
% ylabel('Your Value');
% title('Sparse Frame Points per File');
% 
% % Optional: set x ticks to file indices with labels
% xticks(1:num_files);
% xticklabels(arrayfun(@(x) sprintf('File %d', x), 1:num_files, 'UniformOutput', false));
% % xlim([1, num_files + 0.1]);
% grid on;

