close all; clear all; clc;


libpath = 'U:\Direction_Finding\20251107_Marana_915MHz_CalibrationLibrary_360AZ_66_-6_EL_step3\data';
testpath = 'U:\Direction_Finding\20251112_Marana_915MHz_CalibrationLibrary_360AZ_66_-6_EL_step3_incomplete\data';

az = 66;
el = 48;
f = 0.915;
AZ = 180:-3:-180;
EL = 66:-3:39;

azstep = find(AZ == az);
elstep = find(EL == el);

fileidx = azstep + length(AZ)*(elstep - 1);
 I = 962;


%load lib file
filename = [libpath '\' num2str(fileidx, '%04d') '.BIN'];
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
end
phase_lib = [phase_1;phase_2;phase_3;phase_4;phase_5;phase_6];
a_lib = [a1;a2;a3;a4;a5;a6];


            figure(30)
        
            set(gcf, 'Position', [100, 100, 1400, 700]);
        
        
            subplot(1, 2, 1)
            plot(1:1024, phase_lib.');
            xlim([1 1024]); ylim([-180 180]);
            xlabel('Frame'); ylabel('Phase Difference [deg]'); 
            title(['Phase']);
            legend('A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'location', 'best')
        
            subplot(1, 2, 2)
            plot(1:1024, a_lib.');
            xlim([1 1024]); ylim([80 120]);
            xlabel('Frame'); ylabel('Magnitude [dB]'); 
            title(['Magnitude']);
            legend('A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'location', 'best')
            sgtitle(['Library at AZ = ' num2str(az) ', EL = ' num2str(el)])



filename = [testpath '\' num2str(fileidx, '%04d') '.BIN'];
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
end
phase_test = [phase_1;phase_2;phase_3;phase_4;phase_5;phase_6];
a_test = [a1;a2;a3;a4;a5;a6];

            figure(31)
        
            set(gcf, 'Position', [100, 100, 1400, 700]);
        
        
            subplot(1, 2, 1)
            plot(1:1024, phase_test.');
            xlim([1 1024]); ylim([-180 180]);
            xlabel('Frame'); ylabel('Phase Difference [deg]'); 
            title(['Phase']);
            legend('A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'location', 'best')

            subplot(1, 2, 2)
            plot(1:1024, a_test.');
            xlim([1 1024]); ylim([80 120]);
            xlabel('Frame'); ylabel('Magnitude [dB]'); 
            title(['Magnitude']);
            legend('A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'location', 'best')
            sgtitle(['Test at AZ = ' num2str(az) ', EL = ' num2str(el)])






            figure(32)

            set(gcf, 'Position', [100, 100, 1400, 700]);
       
        
            subplot(1, 2, 1)
            plot(1:1024, (mod((phase_test - phase_lib)+180, 360)-180).');
            xlim([1 1024]); 
            xlabel('Frame'); ylabel('Phase Difference [deg] (Test - Lib)'); 
            title(['Phase']);
            legend('A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'location', 'best')

            subplot(1, 2, 2)
            plot(1:1024, (a_test - a_lib).');
            xlim([1 1024]);
            xlabel('Frame'); ylabel('Magnitude [dB] (Test - Lib)'); 
            title(['Magnitude']);
            legend('A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'location', 'best')
            sgtitle({
                ['Library at AZ = ' num2str(az) ', EL = ' num2str(el)]
                'Frame vs Frame Difference In Test Data vs Library'
            });
