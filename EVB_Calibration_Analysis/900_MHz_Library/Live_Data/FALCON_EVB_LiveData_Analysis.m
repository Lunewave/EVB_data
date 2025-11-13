close all; clear all; clc;

save_figs = 1;
test_location = 'Incomplete Library';
noise_level_test = 45;
testpath = 'U:\Direction_Finding\20251112_Marana_915MHz_CalibrationLibrary_360AZ_66_-6_EL_step3_incomplete\data';

%%%%%%%%%%% FIXED PARAMETERS %%%%%%%%%%%%%
AZ_start = -180; AZ_end = 180; AZ_step = 3;
EL_start = 66; EL_end = -6; EL_step = -3;
AZ_data = AZ_start:AZ_step:AZ_end;
AZ_steps = length(AZ_data);
EL_data = EL_start:EL_step:EL_end;
EL_steps = length(EL_data);
numpeaks2check = 1; %# of peaks to check in each dimension of angle interpolation
%%%%%%%%%%%% LIBRARY %%%%%%%%%%%%%%%%%%%%%
frequency = 915; %MHz
lib_location = 'Calibration Library';
noise_level_cal = 45;
libpath = 'U:\Direction_Finding\20251107_Marana_915MHz_CalibrationLibrary_360AZ_66_-6_EL_step3\data';
offset = 0;
lib_cache = fullfile(libpath, [num2str(frequency/1000) 'GHz_cached_library_data.mat']);
if isfile(lib_cache)
    load(lib_cache, 'Lib_Mag', 'Lib_Phase', 'Lib_Complex');
else
    [Lib_Mag, Lib_Phase, Lib_Complex] = Load_FALCON_EVB_Data_900MHz(libpath, AZ_steps, EL_steps, offset);
    save(lib_cache, 'Lib_Mag', 'Lib_Phase', 'Lib_Complex');
end
%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%
offset = 0;
data_freq = 0.915; %Frequency of test data signal in GHz
test_cache = fullfile(testpath, [num2str(data_freq) 'GHz_cached_test_data.mat']);
if isfile(test_cache)
    load(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex', 'num_files', 'numgoodframes');
else
    [Test_Mag, Test_Phase, Test_Complex, num_files, numgoodframes] = Load_FALCON_EVB_LiveData_900MHz(testpath, offset, data_freq);
    save(test_cache, 'Test_Mag', 'Test_Phase', 'Test_Complex', 'num_files', 'numgoodframes');
end

%% Angle Finding with Interpolation
%%%%%%%%%%%%  DF code  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AZ_table=AZ_data;
EL_table=EL_data;
length_el=length(EL_table);
length_az=length(AZ_table);
AF_results = zeros(num_files, 2);
AF_ITP_results = zeros(num_files, 2);
for ifile= 1:num_files   %-140:2:140;                      % object AZ angle input   [-140:140]
        ifile
        %Test vector
        test=squeeze(Test_Complex(:, ifile)); 

        if isnan(test)
            AF_results(ifile,:) = [NaN NaN];
            AF_ITP_results(ifile,:) = [NaN NaN];
        else
            weighting=[1 1 1 1 1 1];
            for az=1:AZ_steps
                for el=1:EL_steps
                    r=corrcoef(test.*weighting',squeeze(Lib_Complex(:,az,el)).*weighting');
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

% 
% figure(1000)
% set(gcf, 'Position',  [200, 200, 800, 600]);
% if EL_steps == 1
%     plot(AZ_data, mag2db(abs(coe(:,1:end).')))
%     xlabel('AZ');ylabel('Correlation Mag')
% elseif AZ_steps == 1
%     plot(EL_data, mag2db(abs(coe(:,1:end).')))
%     xlabel('EL');ylabel('Correlation Mag')
%     ylim([-2 -0.25])
% else
%     imagesc(AZ_data,EL_data,mag2db(abs(coe(:,1:end).')))
%     caxis([-8 0]);
%     colorbar;
%     xlabel('AZ');ylabel('EL');
% end
% title(['Frame:  ' num2str(ifile)]);
%% Plot Figures
figure(1)
sgtitle(sprintf('Azimuth Results'));
subplot(1, 2, 1)
plot(AF_results(:, 1), 'o-')
title('No Interpolation')
xlabel('File')
ylabel('Azimuth Angle')
grid on
legend('Test Data', 'Ground Truth', 'location', 'best')
subplot(1, 2, 2)
plot(AF_ITP_results(:, 1), 'o-')
title('Interpolation')
xlabel('File')
ylabel('Azimuth Angle')
grid on
legend('Test Data', 'Ground Truth', 'location', 'best')
set(gcf, 'Position', [100, 100, 1400, 700]);


figure(2)
sgtitle(sprintf('Elevation Results'));
subplot(1, 2, 1)
plot(AF_results(:, 2), 'o-')
title('No Interpolation')
xlabel('File')
ylabel('Elevation Angle')
grid on
legend('Test Data', 'Ground Truth', 'location', 'best')
subplot(1, 2, 2)
plot(AF_ITP_results(:, 2), 'o-')
title('Interpolation')
xlabel('File')
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



figure(6)
plot(Test_Mag', 'o-')
ylabel('Magnitude (dB)')
xlabel('File')
title('Signal Strength vs File')
legend('Antenna 1', 'Antenna 2', 'Antenna 3', 'Antenna 4', 'Antenna 5', 'Antenna 6', 'location', 'best')
set(gcf, 'Position', [100, 100, 1400, 700]);





%go from angle to indices
AF_idxs(:, 1) = (AF_results(:, 1) - AZ_start)/(AZ_step) + 1;
AF_idxs(:, 2) = (AF_results(:, 2) - EL_start)/(EL_step) + 1;
[m, n, p] = size(Lib_Mag);
cols   = AF_idxs(:,1);
slices = AF_idxs(:,2);
libmag = zeros(m, num_files);
ax_string = strings(1, num_files);

for k = 1:numel(cols)
    if isnan(cols(k)) || isnan(slices(k))
        % If either index is NaN, store NaN in libmag
        libmag(:, k) = NaN;
        ax_string{k} = '(NaN, NaN)';
    else
        % Otherwise, index into Lib_Mag as normal
        libmag(:, k) = Lib_Mag(:, cols(k), slices(k));
        ax_string{k} = ['(' num2str(AF_results(k, 1)) ',' num2str(AF_results(k, 2)) ')'];
    end
end









figure(101)
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
plot(Test_Mag', 'o-')
grid on
ylabel('Magnitude (dB)')
xlabel('File')
xlim([1 length(AF_results)])
title('Test Signal Strength vs File')
legend('Antenna 1', 'Antenna 2', 'Antenna 3', 'Antenna 4', 'Antenna 5', 'Antenna 6', 'location', 'best')
nexttile
plot(libmag', 'o-')
grid on
ylabel('Magnitude (dB)')
xlabel('File')
xlim([1 length(AF_results)])
title('Library Signal Strength vs File')
legend('Antenna 1', 'Antenna 2', 'Antenna 3', 'Antenna 4', 'Antenna 5', 'Antenna 6', 'location', 'best')
set(gcf, 'Position', [100, 100, 1400, 700]);

figure(102)
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
hold on
colors = lines(6);
markers = {'-o', '-d'};
for i = 1:6
    plot(Test_Mag(i, :), markers{1}, 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    plot(libmag(i, :), markers{2}, 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
end
grid on
ylabel('Magnitude (dB)')
xlabel('Angle (AZ, EL)')
xlim([1 length(AF_results)])
legend({'Test A1','Lib A1','Test A2','Lib A2','Test A3','Lib A3', ...
        'Test A4','Lib A4','Test A5','Lib A5','Test A6','Lib A6'}, ...
        'Location','best','NumColumns',6)
xticks(1:num_files)
xticklabels(ax_string)
xtickangle(45)   % tilt for readability
set(gcf, 'Position', [100, 100, 1400, 700]);

testdif = zeros(5, num_files);
libdif = zeros(5, num_files);

for k = 1:5
    testdif(k, :) = Test_Mag(k, :) - Test_Mag(6, :);
    libdif(k, :) = libmag(k, :) - libmag(6, :);
end
test_average = mean(testdif, 2, 'omitnan');
lib_average = mean(libdif, 2, 'omitnan');


diff = testdif - libdif;


figure(103)
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Test signal differences ---
nexttile
plot(testdif', 'o-')
grid on
ylabel('Magnitude (dB)')
xlabel('File')
xlim([1 length(AF_results)])
title('Test Signal Magnitude Difference vs File')

% Build legend as cell array
test_legend = {
    sprintf('A1 - A6, Mean: %.2f', test_average(1)), 
    sprintf('A2 - A6, Mean: %.2f', test_average(2)), 
    sprintf('A3 - A6, Mean: %.2f', test_average(3)), 
    sprintf('A4 - A6, Mean: %.2f', test_average(4)), 
    sprintf('A5 - A6, Mean: %.2f', test_average(5))
};
legend(test_legend, 'Location', 'best')

% --- Library signal differences ---
nexttile
plot(libdif', 'o-')
grid on
ylabel('Magnitude (dB)')
xlabel('File')
xlim([1 length(AF_results)])
title('Library Signal Magnitude Difference vs File')

% Build legend as cell array
lib_legend = {
    sprintf('A1 - A6, Mean: %.2f', lib_average(1)), 
    sprintf('A2 - A6, Mean: %.2f', lib_average(2)), 
    sprintf('A3 - A6, Mean: %.2f', lib_average(3)), 
    sprintf('A4 - A6, Mean: %.2f', lib_average(4)), 
    sprintf('A5 - A6, Mean: %.2f', lib_average(5))
};
legend(lib_legend, 'Location', 'best')

set(gcf, 'Position', [100, 100, 1400, 700]);



figure(104)
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
plot(diff', 'o-')
grid on
ylabel('Magnitude (dB)')
xlabel('Angle (AZ, EL)')
xlim([1 length(AF_results)])
legend('Test(A1 - A6) - Lib(A1 - A6)', 'Test(A2 - A6) - Lib(A2 - A6)', 'Test(A3 - A6) - Lib(A3 - A6)', 'Test(A4 - A6) - Lib(A4 - A6)', 'Test(A5 - A6) - Lib(A5 - A6)')
xticks(1:num_files)
xticklabels(ax_string)
xtickangle(45)   % tilt for readability
set(gcf, 'Position', [100, 100, 1400, 700]);

num_ant = size(Test_Mag, 1);   % number of antennas
num_files = size(Test_Mag, 2);

% Initialize: residuals(i,j,k) = difference between test vs lib
residuals = zeros(num_ant, num_ant, num_files);

for i = 1:num_ant
    for j = 1:num_ant
        if i ~= j
            % Relative diff between test pair and lib pair
            test_diff = Test_Mag(i,:) - Test_Mag(j,:);
            lib_diff  = libmag(i,:)  - libmag(j,:);
            residuals(i,j,:) = test_diff - lib_diff;
        end
    end
end
rms_residuals = squeeze(mean(residuals, 3, 'omitnan')); % antenna × antenna matrix
figure(201)
h = heatmap(1:num_ant, 1:num_ant, rms_residuals, ...
    'Colormap', parula, 'ColorbarVisible','on');
h.XLabel = 'Antenna j';
h.YLabel = 'Antenna i';
h.Title = 'Signed Average Difference: (Test(i)-Test(j)) - (Lib(i)-Lib(j))';
h.ColorLimits = [-5 5];






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
    saveas(figure(4), fullfile(newFolderPath, 'Good_Frames.jpeg'));
    saveas(figure(6), fullfile(newFolderPath, 'Signal_Strength.jpeg'));
    saveas(figure(101), fullfile(newFolderPath, 'Signal_Strength_vs_Lib.jpeg'));
    saveas(figure(102), fullfile(newFolderPath, 'Signal_Strength_vs_Lib_1fig.jpeg'));
    saveas(figure(103), fullfile(newFolderPath, 'libdif_testdif.jpeg'));
    saveas(figure(104), fullfile(newFolderPath, 'A6_residuals.jpeg'));
    saveas(figure(201), fullfile(newFolderPath, 'residuals.jpeg'));




    parentFolder = fileparts(newFolderPath); % get parent folder
    save(fullfile(parentFolder, 'AF_ITP_results.mat'), 'AF_ITP_results');
end





test_lib_p = zeros(6, 121, 10);
test_lib_m = zeros(6, 121, 10);
for i =1:10
    test_lib_p(:, :, i) = Test_Phase(:, (1:121)+121*(i-1));
    test_lib_m(:, :, i) = Test_Mag(:, (1:121)+121*(i-1));
end


for i = 1:6
    figure(i)
    t = tiledlayout(1, 2, 'Padding', 'tight', 'TileSpacing', 'tight'); % Tight spacing

    % Compute global min/max for this antenna
    data1 = abs(squeeze(Lib_Phase(i, :, 1:10))).';
    data2 = abs(squeeze(test_lib_p(i, :, :))).';
    cmin = min([data1(:); data2(:)]);
    cmax = max([data1(:); data2(:)]);

    % Calibration Library
    nexttile
    imagesc(-180:3:180, 66:-3:39, data1);
    caxis([cmin-5 cmax+5]);
    cb1 = colorbar;
    cb1.Label.String = 'Phase [deg]';
    xlabel('AZ'); ylabel('EL'); 
    title('Calibration Library');

    % Test Library
    nexttile
    imagesc(-180:3:180, 66:-3:39, data2);
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
    data1 = abs(squeeze(Lib_Mag(i, :, 1:10))).';
    data2 = abs(squeeze(test_lib_m(i, :, :))).';
    cmin = min([data1(:); data2(:)]);
    cmax = max([data1(:); data2(:)]);

    % Calibration Library
    nexttile
    imagesc(-180:3:180, 66:-3:39, data1);
    caxis([cmin-5 cmax+5]);
    cb1 = colorbar;
    cb1.Label.String = 'Magnitude [dB]';
    xlabel('AZ'); ylabel('EL'); 
    title('Calibration Library');

    % Test Library
    nexttile
    imagesc(-180:3:180, 66:-3:39, data2);
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
    data1 = abs(squeeze(Lib_Phase(i, :, 1:10))).';
    data2 = abs(squeeze(test_lib_p(i, :, :))).';
    data = abs(data2 - data1);
    cmin = min(data(:));
    cmax = max(data(:));

    % Calibration Library
    nexttile
    imagesc(-180:3:180, 66:-3:39, data);
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
    data1 = abs(squeeze(Lib_Mag(i, :, 1:10))).';
    data2 = abs(squeeze(test_lib_m(i, :, :))).';
    data = abs(data2 - data1);
    cmin = min(data(:));
    cmax = max(data(:));

    % Calibration Library
    nexttile
    imagesc(-180:3:180, 66:-3:39, data);
    caxis([cmin cmax+1]);
    cb1 = colorbar;
    cb1.Label.String = 'Magnitude [dB]';
    xlabel('AZ'); ylabel('EL'); 
    title(['Antenna ' num2str(i) ' Test Library - Calibration Library (abs(Magnitude))']);

    set(gcf, 'Position', [100, 100, 1400, 700]);


end

DF = zeros(2, 121, 10);
for i =1:10
    DF(1, :, i) = AF_ITP_results((1:121)+121*(i-1), 1);
    DF(2, :, i) = AF_ITP_results((1:121)+121*(i-1), 2);
end

for i = 1:10
    DF(1, :, i) = DF(1, :, i) - (-180:3:180);
    DF(2, :, i) = DF(2, :, i) - 66 + 3*(i-1);
end

DF(1, :, :)=mod(DF(1, :, :)+180,360)-180;
DF(2, :, :)=mod(DF(2, :, :)+180,360)-180;


figure(1001)
t = tiledlayout(1, 1, 'Padding', 'tight', 'TileSpacing', 'tight'); % Tight spacing

% Compute global min/max for this antenna
data1 = abs(squeeze(DF(1, :, :))).';
data2 = abs(squeeze(DF(2, :, :))).';
% cmin = min([data1(:); data2(:)]);
% cmax = max([data1(:); data2(:)]);

% Calibration Library
nexttile
imagesc(-180:3:180, 66:-3:39, data1);
% caxis([cmin-5 cmax+5]);
cb1 = colorbar;
cb1.Label.String = 'Azimuth Error [deg]';
xlabel('AZ'); ylabel('EL'); 
title('Interpolated Azimuth Error');

% % Test Library
% nexttile
% imagesc(-180:3:180, 66:-3:39, data2);
% % caxis([cmin-5 cmax+5]);
% cb2 = colorbar;
% cb2.Label.String = 'Elevation Error [deg]';
% xlabel('AZ'); ylabel('EL'); 
% title('Interpolated Elevation Error');

set(gcf, 'Position', [100, 100, 1400, 700]);
sgtitle(t, ['Angle Finding Results'])


figure(1002)
t = tiledlayout(1, 1, 'Padding', 'tight', 'TileSpacing', 'tight'); % Tight spacing

% Compute global min/max for this antenna
data1 = abs(squeeze(DF(1, :, :))).';
data2 = abs(squeeze(DF(2, :, :))).';
% cmin = min([data1(:); data2(:)]);
% cmax = max([data1(:); data2(:)]);

% % Calibration Library
% nexttile
% imagesc(-180:3:180, 66:-3:39, data1);
% % caxis([cmin-5 cmax+5]);
% cb1 = colorbar;
% cb1.Label.String = 'Azimuth Error [deg]';
% xlabel('AZ'); ylabel('EL'); 
% title('Interpolated Azimuth Error');

% Test Library
nexttile
imagesc(-180:3:180, 66:-3:39, data2);
% caxis([cmin-5 cmax+5]);
cb2 = colorbar;
cb2.Label.String = 'Elevation Error [deg]';
xlabel('AZ'); ylabel('EL'); 
title('Interpolated Elevation Error');

set(gcf, 'Position', [100, 100, 1400, 700]);
sgtitle(t, ['Angle Finding Results'])