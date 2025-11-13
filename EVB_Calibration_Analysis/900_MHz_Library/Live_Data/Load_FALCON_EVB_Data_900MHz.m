function [magnitude, phase, complex_values] = Load_FALCON_EVB_Data(path,numAZ, numEL, offset, data_freq)
% Load_FALCON_EVB_Data loads and processes binary EVB data into magnitude, phase, and complex values.
%
% Inputs:
%   path   - Folder path to .BIN files
%   numAZ  - Number of azimuth steps
%   numEL  - Number of elevation steps
%
% Outputs:
%   magnitude      - [6 x numAZ x numEL] matrix of dB values
%   phase          - [6 x numAZ x numEL] matrix of phase in degrees
%   complex_values - [6 x numAZ x numEL] matrix of complex values (linear scale)

    magnitude = zeros(6, numAZ, numEL);
    phase = zeros(6, numAZ, numEL);
    complex_values = zeros(6, numAZ, numEL);
    fc_MHz = 900;
    fre_sample=2.94912e9/12
    bw_MHz = fre_sample/1e6
    fre=[0:1023]/1024*fre_sample;
    a = fre/1e9+(fc_MHz/1e3 - bw_MHz/(2*1e3));
    n = length(fre);
    a = a(end:-1:1);
    a = [a(n/2+1:end),a(1:n/2)];

    [~, I] = min(abs(a-data_freq));
    I = I+1;
    for i = 1:numEL
        for k = 1:numAZ
            file_number = (i-1)*numAZ + k-1 + offset
            filename = [path '\' num2str(file_number,'%04d') '.BIN']
            fileID = fopen(filename, 'r', 'ieee-le');
            C = fread(fileID, Inf, 'int16');fclose(fileID);
            C0 = reshape(C,[2,length(C)/2]).';
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
            
            for frame_ind=1:1024
                % freA=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),1));
                freB=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),2)); %CH2
                freC=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),3)); %CH3
                freD=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),4)); %CH4
                freE=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),5)); %CH5
                freF=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),6)); %CH6
                freG=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),7)); %CH7
                % freH=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),8)); %CH8


                % [~,I]=max(abs(freD));
                % I=234; %2.456 GHz
                % Phase Difference
                phase_5(frame_ind)=angle(freB(I)/freE(I))/pi*180; %EVB 2 vs 5 i.e phase difference of antenna 5 relative to antenna 1
                phase_4(frame_ind)=angle(freC(I)/freE(I))/pi*180; %EVB 3 vs 5 i.e phase difference of antenna 4 relative to antenna 1
                phase_1(frame_ind)=angle(freE(I)/freE(I))/pi*180; %EVB 5 vs 5 i.e phase difference of antenna 1 relative to antenna 1
                phase_2(frame_ind)=angle(freF(I)/freE(I))/pi*180; %EVB 6 vs 5 i.e phase difference of antenna 2 relative to antenna 1
                phase_3(frame_ind)=angle(freG(I)/freE(I))/pi*180; %EVB 7 vs 5 i.e phase difference of antenna 3 relative to antenna 1
                phase_6(frame_ind)=angle(freD(I)/freE(I))/pi*180; %EVB 4 vs 5 i.e phase difference of antenna 6 relative to antenna 1


                % phase_5(frame_ind)=angle(freB(I)/freD(I))/pi*180;
                % phase_4(frame_ind)=angle(freC(I)/freD(I))/pi*180;
                % phase_1(frame_ind)=angle(freE(I)/freH(I))/pi*180;
                % phase_2(frame_ind)=angle(freF(I)/freH(I))/pi*180;
                % phase_3(frame_ind)=angle(freG(I)/freH(I))/pi*180;
                % phase_6(frame_ind)=0;





                % Magnitude
                c5(frame_ind) = mag2db(abs(freB(I)));
                c4(frame_ind) = mag2db(abs(freC(I)));
                c1(frame_ind) = mag2db(abs(freE(I)));
                c2(frame_ind) = mag2db(abs(freF(I)));
                c3(frame_ind) = mag2db(abs(freG(I)));
                c6(frame_ind) = mag2db(abs(freD(I)));

            end
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
end

