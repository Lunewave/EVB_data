close all; clear all; clc;
AZ_start = -180; AZ_end = 180; AZ_step = 3;
EL_start = 66; EL_end = 0; EL_step = -3;


AZ_data = AZ_start:AZ_step:AZ_end;
numAZ = length(AZ_data);
EL_data = EL_start:EL_step:EL_end;
numEL = length(EL_data);
offset = 0;

path = 'U:\Falcon_Project\20250625_MaranaTest_AZ360_EL66_Step3_withLens_withEVB_2.456GHz_CalibrationLibrary';

magnitude = zeros(6, numAZ, numEL);
phase = zeros(6, numAZ, numEL);
complex_values = zeros(6, numAZ, numEL);
percent_table = zeros(numAZ, numEL);

for i = 1:numEL
    for k = 1:numAZ
        num2str((i-1)*numAZ + k-1,'%04d')
        fileID = fopen([path '\' num2str((i-1)*numAZ + k-1 + offset,'%04d') '.BIN'], 'r', 'ieee-le');
        C = fread(fileID, Inf, 'int16');fclose(fileID);
        C0 = reshape(C,[8,length(C)/8]).';
        
        % C1= C0 (1:8:end,:).'; C_all(:,1)=C1(:);
        C2= C0 (2:8:end,:).'; C_all(:,2)=C2(:);
        C3= C0 (3:8:end,:).'; C_all(:,3)=C3(:);
        C4= C0 (4:8:end,:).'; C_all(:,4)=C4(:);
        C5= C0 (5:8:end,:).'; C_all(:,5)=C5(:);
        C6= C0 (6:8:end,:).'; C_all(:,6)=C6(:);
        C7= C0 (7:8:end,:).'; C_all(:,7)=C7(:);
        C8= C0 (8:8:end,:).'; C_all(:,8)=C8(:); 
        
        C1_cmplex=C_all(1:2:end,:)+1i*C_all(2:2:end,:);

        
        frame_valid = zeros(1, 512);  % 0 or 1 for each frame
        threshold_dB = 65;

        for frame_ind=1:512
            % freA=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),1));
            freB=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),2)); %CH2
            freC=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),3)); %CH3
            freD=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),4)); %CH4
            freE=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),5)); %CH5
            freF=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),6)); %CH6
            freG=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),7)); %CH7
            freH=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),8)); %CH8
            % Get dB values at index I = 844
            [~,I]=max(abs(freD));
            I=844; %2.356 GHz
            dB_vals = [
                mag2db(abs(freB(I))),  % Antenna 5
                mag2db(abs(freC(I))),  % Antenna 4
                mag2db(abs(freD(I))),  % Antenna 6
                mag2db(abs(freE(I))),  % Antenna 1 (ref)
                mag2db(abs(freF(I))),  % Antenna 2
                mag2db(abs(freG(I)))   % Antenna 3
            ];
            
            % Check if all antennas exceed the threshold
            if all(dB_vals > threshold_dB)
                frame_valid(frame_ind) = 1;
            end

            
           

            % Phase Difference
            phase_5(frame_ind)=angle(freB(I)/freE(I))/pi*180; %EVB 2 vs 5 i.e phase difference of antenna 5 relative to antenna 1
            phase_4(frame_ind)=angle(freC(I)/freE(I))/pi*180; %EVB 3 vs 5 i.e phase difference of antenna 4 relative to antenna 1
            phase_1(frame_ind)=angle(freE(I)/freE(I))/pi*180; %EVB 5 vs 5 i.e phase difference of antenna 1 relative to antenna 1
            phase_2(frame_ind)=angle(freF(I)/freE(I))/pi*180; %EVB 6 vs 5 i.e phase difference of antenna 2 relative to antenna 1
            phase_3(frame_ind)=angle(freG(I)/freE(I))/pi*180; %EVB 7 vs 5 i.e phase difference of antenna 3 relative to antenna 1
            phase_6(frame_ind)=angle(freD(I)/freE(I))/pi*180; %EVB 4 vs 5 i.e phase difference of antenna 6 relative to antenna 1

            % Magnitude
            c5(frame_ind) = mag2db(abs(freB(I)));
            c4(frame_ind) = mag2db(abs(freC(I)));
            c1(frame_ind) = mag2db(abs(freE(I)));
            c2(frame_ind) = mag2db(abs(freF(I)));
            c3(frame_ind) = mag2db(abs(freG(I)));
            % if mag2db(abs(freD(I)))>=mag2db(abs(freH(I)))
            c6(frame_ind) = mag2db(abs(freD(I)));
            % else
            %     c6(frame_ind) = mag2db(abs(freH(I)));
            % end
        end
        percent_valid = 100 * sum(frame_valid) / length(frame_valid);
        fprintf('Percentage of valid frames: %.2f%%\n', percent_valid);
        percent_table(k, i) = percent_valid;





        P_diff1=(mod(phase_1+180,360)-180).';
        P_diff2=(mod(phase_2+180,360)-180).';
        P_diff3=(mod(phase_3+180,360)-180).';
        P_diff4=(mod(phase_4+180,360)-180).';
        P_diff5=(mod(phase_5+180,360)-180).';
        P_diff6=(mod(phase_6+180,360)-180).';

        phase(:,k, i)=[median(P_diff1) median(P_diff2) median(P_diff3) median(P_diff4) median(P_diff5) median(P_diff6)];
        magnitude(:, k, i) = [median(c1) median(c2) median(c3) median(c4) median(c5) median(c6)];
        for antenna = 1:6
            complex_values(antenna, k, i) = db2mag(magnitude(antenna, k, i))*exp(j.*phase(antenna, k, i)/180*pi);
        end
    end
end

save('processed_data.mat', 'phase', 'magnitude', 'complex_values');


figure;
imagesc(AZ_data, EL_data, percent_table');
set(gca, 'YDir', 'normal');  % Flip Y so EL=0 is at bottom
colorbar;
xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');
title('Percent of Valid Frames (All Antennas > 65 dB)');
clim([0 100])


