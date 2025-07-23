close all; clear all; clc;

save_figs = 1;
frequency = 2456; %MHz
test_location = 'Drone Test';
noise_level_test = 45;
testpath = 'U:\Falcon_Project\20250711_MaranaTest_AZ360_EL0_Step5_withLens_withEVB_2.456GHz_DroneTest_r-10_h-2';
num_frames = 1024;
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
data_freq = 2.357;
fre_sample=2.94912e9/12;
fre=[0:1023]/1024*fre_sample;
a = fre/1e9+2.277;
n = length(fre);
a = [a(n/2+1:end),a(1:n/2)];
a = a-data_freq;
[~, I] = min(abs(a));
bin_files = dir(fullfile(testpath, '*.BIN'));
num_files = numel(bin_files)-1;

for file = 1:num_files
num2str(file,'%04d')
fileID = fopen([testpath '\' num2str(file,'%04d') '.BIN'], 'r', 'ieee-le');
C = fread(fileID, Inf, 'int16');fclose(fileID);
C0 = reshape(C,[8,length(C)/8]).';
L_C0=length(C0);

C1= C0 (1:2:L_C0/4,:).'; C_all(:,1)=C1(:);
C2= C0 (2:2:L_C0/4,:).'; C_all(:,2)=C2(:);
C3= C0 (L_C0/4+1:2:L_C0/4*2,:).'; C_all(:,3)=C3(:);
C4= C0 (L_C0/4+2:2:L_C0/4*2,:).'; C_all(:,4)=C4(:);
C5= C0 (L_C0/2+1:2:L_C0/4*3,:).'; C_all(:,5)=C5(:);
C6= C0 (L_C0/2+2:2:L_C0/4*3,:).'; C_all(:,6)=C6(:);
C7= C0 (L_C0/4*3+1:2:L_C0,:).'; C_all(:,7)=C7(:);
C8= C0 (L_C0/4*3+2:2:L_C0,:).'; C_all(:,8)=C8(:); 

C1_cmplex=C_all(1:2:end,:)+1i*C_all(2:2:end,:);

% figure(1000)
for frame_ind=1:num_frames
    % freA=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),1));
    freB=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),2)); %CH2
    freC=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),3)); %CH3
    freD=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),4)); %CH4
    freE=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),5)); %CH5
    freF=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),6)); %CH6
    freG=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),7)); %CH7
    % freH=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),8)); %CH8
    
    
    % Magnitude
    a5 = mag2db(abs(freB(I)));
    a4 = mag2db(abs(freC(I)));
    a1 = mag2db(abs(freE(I)));
    a2 = mag2db(abs(freF(I)));
    a3 = mag2db(abs(freG(I)));
    a6 = mag2db(abs(freD(I)));


    % Phase Difference
    phase_5=angle(freB(I)/freE(I))/pi*180; %EVB 2 vs 5 i.e phase difference of antenna 5 relative to antenna 1
    phase_4=angle(freC(I)/freE(I))/pi*180; %EVB 3 vs 5 i.e phase difference of antenna 4 relative to antenna 1
    phase_1=angle(freE(I)/freE(I))/pi*180; %EVB 5 vs 5 i.e phase difference of antenna 1 relative to antenna 1
    phase_2=angle(freF(I)/freE(I))/pi*180; %EVB 6 vs 5 i.e phase difference of antenna 2 relative to antenna 1
    phase_3=angle(freG(I)/freE(I))/pi*180; %EVB 7 vs 5 i.e phase difference of antenna 3 relative to antenna 1
    phase_6=angle(freD(I)/freE(I))/pi*180; %EVB 4 vs 5 i.e phase difference of antenna 6 relative to antenna 1


    signal_threshold = 60;
    
    if a1 < signal_threshold || ...
       a2 < signal_threshold || ...
       a3 < signal_threshold || ...
       a4 < signal_threshold || ...
       a5 < signal_threshold || ...
       a6 < signal_threshold
    
        a1 = NaN;
        a2 = NaN;
        a3 = NaN;
        a4 = NaN;
        a5 = NaN;
        a6 = NaN;
    
        phase_1 = NaN;
        phase_2 = NaN;
        phase_3 = NaN;
        phase_4 = NaN;
        phase_5 = NaN;
        phase_6 = NaN;
    end
    P_diff1=(mod(phase_1+180,360)-180).';
    P_diff2=(mod(phase_2+180,360)-180).';
    P_diff3=(mod(phase_3+180,360)-180).';
    P_diff4=(mod(phase_4+180,360)-180).';
    P_diff5=(mod(phase_5+180,360)-180).';
    P_diff6=(mod(phase_6+180,360)-180).';
    magnitude = [a1 a2 a3 a4 a5 a6];
    phase = [P_diff1 P_diff2 P_diff3 P_diff4 P_diff5 P_diff6];

    for antenna = 1:6
        Test_Complex(antenna, frame_ind) = db2mag(magnitude(antenna))*exp(j.*phase(antenna)/180*pi);
        % subplot(2, 3, antenna)
        % scatter(frame_ind, magnitude(antenna))
        % hold on
    end
end

%% Angle Finding with Interpolation
%%%%%%%%%%%%  DF code  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AZ_table=AZ_data;
EL_table=[EL_start:EL_step:EL_end];
length_el=length(EL_table);
length_az=length(AZ_table);
AF_results = zeros(num_frames, 2);
AF_ITP_results = zeros(num_frames, 2);
for ifile= 1:num_frames   %-140:2:140;                      % object AZ angle input   [-140:140]

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
Median_Results(file, :) = [median(AF_results(:, 1), 'omitnan') median(AF_results(:, 2), 'omitnan')];
Median_ITP_Results(file, :) = [median(AF_ITP_results(:, 1), 'omitnan') median(AF_ITP_results(:, 2), 'omitnan')];

end



figure(1)
subplot(1, 2, 1)
scatter(1:num_files, Median_ITP_Results(:, 1))
hold on
scatter(1:num_files, Median_Results(:, 1))
hold off
xlabel('Frame')
ylabel('Azimuth DF Result')
xlim([0 num_files])
legend('Interpolated', 'Not Interpolated', 'location', 'best')

subplot(1, 2, 2)
scatter(1:num_files, Median_ITP_Results(:, 2))
hold on
scatter(1:num_files, Median_Results(:, 2))
hold off
xlabel('Frame')
ylabel('Elevation DF Result')
xlim([0 num_files])
legend('Interpolated', 'Not Interpolated', 'location', 'best')
set(gcf, 'Position', [100, 100, 1400, 700]);











