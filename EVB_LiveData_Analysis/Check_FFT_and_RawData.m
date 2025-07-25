%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all;close all;clc;
path = 'U:\Falcon_Project\20250711_MaranaTest_AZ360_EL0_Step5_withLens_withEVB_2.456GHz_DroneTest_r-10_h-2\';
frame = 36;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen([path num2str(frame,'%04d') '.BIN'], 'r', 'ieee-le');


% C = fread(fileID, Inf, 'int16');fclose(fileID);
% C0 = reshape(C,[8,length(C)/8]).';
% 
% C1= C0 (1:8:end,:).'; C_all(:,1)=C1(:);
% C2= C0 (2:8:end,:).'; C_all(:,2)=C2(:);
% C3= C0 (3:8:end,:).'; C_all(:,3)=C3(:);
% C4= C0 (4:8:end,:).'; C_all(:,4)=C4(:);
% C5= C0 (5:8:end,:).'; C_all(:,5)=C5(:);
% C6= C0 (6:8:end,:).'; C_all(:,6)=C6(:);
% C7= C0 (7:8:end,:).'; C_all(:,7)=C7(:);
% C8= C0 (8:8:end,:).'; C_all(:,8)=C8(:); 
% 
% C1_cmplex=C_all(1:2:end,:)+1i*C_all(2:2:end,:);

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



figure(1)
set(gcf, 'Position', [100, 100, 1200, 500]);
subplot(1,2,1)
plot(real(C1_cmplex(:,2)));
hold on;
plot(real(C1_cmplex(:,3)));
plot(real(C1_cmplex(:,4)));
plot(real(C1_cmplex(:,5)));
plot(real(C1_cmplex(:,6)));
plot(real(C1_cmplex(:,7)));
xlabel('Sample index');
ylabel('ADC output')
xlim([1 128])
title(['Real Component of ADC Frame ' num2str(frame)]);
legend('Antenna 5','Antenna 4', 'Antenna 6', 'Antenna 1', 'Antenna 2', 'Antenna 3', 'Location', 'best');

subplot(1,2,2)
plot(abs(C1_cmplex(:,2)), 'linewidth', 1.5);
hold on;
plot(abs(C1_cmplex(:,3)), 'linewidth', 1.5);
plot(abs(C1_cmplex(:,4)), 'linewidth', 1.5);
plot(abs(C1_cmplex(:,5)), 'linewidth', 1.5);
plot(abs(C1_cmplex(:,6)), 'linewidth', 1.5);
plot(abs(C1_cmplex(:,7)), 'linewidth', 1.5);
xlabel('Sample index');
ylabel('ADC output')
xlim([1 128])
title(['Absolute Value of ADC Frame ' num2str(frame)]);
legend('Antenna 5','Antenna 4', 'Antenna 6', 'Antenna 1', 'Antenna 2', 'Antenna 3', 'Location', 'best');





fre_sample=2.94912e9/12;
fre=[0:1023]/1024*fre_sample;
for frame_ind=1:512
freA=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),1));
freB=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),2));
freC=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),3));
freD=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),4));

freE=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),5));
freF=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),6));
freG=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),7));
freH=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),8));


n = length(fre);
figure(20)
hold on; grid on;
plot(fre/1e9+2.277, [mag2db(abs(freB(n/2+1:end))); mag2db(abs(freB(1:n/2)))]);
plot(fre/1e9+2.277, [mag2db(abs(freC(n/2+1:end))); mag2db(abs(freC(1:n/2)))]);
plot(fre/1e9+2.277, [mag2db(abs(freD(n/2+1:end))); mag2db(abs(freD(1:n/2)))]);
plot(fre/1e9+2.277, [mag2db(abs(freE(n/2+1:end))); mag2db(abs(freE(1:n/2)))]);
plot(fre/1e9+2.277, [mag2db(abs(freF(n/2+1:end))); mag2db(abs(freF(1:n/2)))]);
plot(fre/1e9+2.277, [mag2db(abs(freG(n/2+1:end))); mag2db(abs(freG(1:n/2)))]);
% plot(fre/1e9+2.277, [mag2db(abs(freH(n/2+1:end))); mag2db(abs(freH(1:n/2)))]);
xlabel('Frequency (GHz)')
ylabel('Magnitude (dB)')
ylim([40 100])
xlim([2.4 2.5])
end


n = length(fre);
figure(2)
hold on; grid on;
plot(fre/1e9+2.277, [mag2db(abs(freB(n/2+1:end))); mag2db(abs(freB(1:n/2)))]);
plot(fre/1e9+2.277, [mag2db(abs(freC(n/2+1:end))); mag2db(abs(freC(1:n/2)))]);
plot(fre/1e9+2.277, [mag2db(abs(freD(n/2+1:end))); mag2db(abs(freD(1:n/2)))]);
plot(fre/1e9+2.277, [mag2db(abs(freE(n/2+1:end))); mag2db(abs(freE(1:n/2)))]);
plot(fre/1e9+2.277, [mag2db(abs(freF(n/2+1:end))); mag2db(abs(freF(1:n/2)))]);
plot(fre/1e9+2.277, [mag2db(abs(freG(n/2+1:end))); mag2db(abs(freG(1:n/2)))]);
xlabel('Frequency (GHz)')
ylabel('Magnitude (dB)')
title(['FFT of ADC Frame ' num2str(frame)]);
legend('Antenna 5','Antenna 4', 'Antenna 6', 'Antenna 1', 'Antenna 2', 'Antenna 3', 'Location', 'best');


