function [magnitude, phase, complex_values, num_files, frames_passed] = Load_FALCON_EVB_LiveData(path,offset, data_freq)
% Load_FALCON_EVB_Data loads and processes binary EVB data into magnitude, phase, and complex values.
%
% Inputs:
%   path   - Folder path to .BIN files

%
% Outputs:
%   magnitude      - [6 x #of bin files] matrix of dB values
%   phase          - [6 x #of bin files] matrix of phase in degrees
%   complex_values - [6 x #of bin files] matrix of complex values
    bin_files = dir(fullfile(path, '*.BIN'));
    num_files = numel(bin_files)-1;


    magnitude = zeros(6, num_files); % 6xnum_files
    phase = zeros(6, num_files);
    complex_values = zeros(6, num_files);
    frames_passed = zeros(1, num_files);
    fre_sample=2.94912e9/12;
    fre=[0:1023]/1024*fre_sample;
    a = fre/1e9+2.277;
    n = length(fre);
    a = [a(n/2+1:end),a(1:n/2)];
    a = a-data_freq;
    [~, I] = min(abs(a));

    for i = 1:num_files %skip first file, bad frame
        num2str(i+offset,'%04d')
        fileID = fopen([path '\' num2str(i+offset,'%04d') '.BIN'], 'r', 'ieee-le');
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

        
        passed = 0;
        for frame_ind=1:1024
            % freA=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),1));
            freB=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),2)); %CH2
            freC=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),3)); %CH3
            freD=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),4)); %CH4
            freE=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),5)); %CH5
            freF=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),6)); %CH6
            freG=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),7)); %CH7
            % freH=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),8)); %CH8
            
            
            % Magnitude
            a5(frame_ind) = mag2db(abs(freB(I)));
            a4(frame_ind) = mag2db(abs(freC(I)));
            a1(frame_ind) = mag2db(abs(freE(I)));
            a2(frame_ind) = mag2db(abs(freF(I)));
            a3(frame_ind) = mag2db(abs(freG(I)));
            a6(frame_ind) = mag2db(abs(freD(I)));


            % Phase Difference
            phase_5(frame_ind)=angle(freB(I)/freE(I))/pi*180; %EVB 2 vs 5 i.e phase difference of antenna 5 relative to antenna 1
            phase_4(frame_ind)=angle(freC(I)/freE(I))/pi*180; %EVB 3 vs 5 i.e phase difference of antenna 4 relative to antenna 1
            phase_1(frame_ind)=angle(freE(I)/freE(I))/pi*180; %EVB 5 vs 5 i.e phase difference of antenna 1 relative to antenna 1
            phase_2(frame_ind)=angle(freF(I)/freE(I))/pi*180; %EVB 6 vs 5 i.e phase difference of antenna 2 relative to antenna 1
            phase_3(frame_ind)=angle(freG(I)/freE(I))/pi*180; %EVB 7 vs 5 i.e phase difference of antenna 3 relative to antenna 1
            phase_6(frame_ind)=angle(freD(I)/freE(I))/pi*180; %EVB 4 vs 5 i.e phase difference of antenna 6 relative to antenna 1


            signal_threshold = 65;
            
            if a1(frame_ind) < signal_threshold || ...
               a2(frame_ind) < signal_threshold || ...
               a3(frame_ind) < signal_threshold || ...
               a4(frame_ind) < signal_threshold || ...
               a5(frame_ind) < signal_threshold || ...
               a6(frame_ind) < signal_threshold
            
                a1(frame_ind) = NaN;
                a2(frame_ind) = NaN;
                a3(frame_ind) = NaN;
                a4(frame_ind) = NaN;
                a5(frame_ind) = NaN;
                a6(frame_ind) = NaN;
            
                phase_1(frame_ind) = NaN;
                phase_2(frame_ind) = NaN;
                phase_3(frame_ind) = NaN;
                phase_4(frame_ind) = NaN;
                phase_5(frame_ind) = NaN;
                phase_6(frame_ind) = NaN;
            else
                passed = passed + 1;
            end

        end
        frames_passed(i) = passed;
        P_diff1=(mod(phase_1+180,360)-180).';
        P_diff2=(mod(phase_2+180,360)-180).';
        P_diff3=(mod(phase_3+180,360)-180).';
        P_diff4=(mod(phase_4+180,360)-180).';
        P_diff5=(mod(phase_5+180,360)-180).';
        P_diff6=(mod(phase_6+180,360)-180).';




        phase(:,i)=[median(P_diff1, 'omitnan') median(P_diff2, 'omitnan') median(P_diff3, 'omitnan') median(P_diff4, 'omitnan') median(P_diff5, 'omitnan') median(P_diff6, 'omitnan')];
        magnitude(:, i) = [median(a1, 'omitnan') median(a2, 'omitnan') median(a3, 'omitnan') median(a4, 'omitnan') median(a5, 'omitnan') median(a6, 'omitnan')];
        for antenna = 1:6
            complex_values(antenna, i) = db2mag(magnitude(antenna, i))*exp(j.*phase(antenna, i)/180*pi);
        end

    end
end
