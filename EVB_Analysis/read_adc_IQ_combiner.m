clear all;close all;clc;
FOV = 6;
path = 'D:\DATA\'; %keep trailing backslash

start_angle = -3;
end_angle = 3;
step = 3;




phase_diff_all = zeros(8, FOV+1);


for ang_ind=start_angle:step:end_angle
fileID = fopen([path num2str((ang_ind+FOV/2)/step,'%04d') '.BIN'], 'r', 'ieee-le');
ang_ind


C = fread(fileID, Inf, 'int16');fclose(fileID);
C0 = reshape(C,[8,length(C)/8]).';

C1= C0 (1:8:end,:).'; C_all(:,1)=C1(:);
C2= C0 (2:8:end,:).'; C_all(:,2)=C2(:);
C3= C0 (3:8:end,:).'; C_all(:,3)=C3(:);
C4= C0 (4:8:end,:).'; C_all(:,4)=C4(:);
C5= C0 (5:8:end,:).'; C_all(:,5)=C5(:);
C6= C0 (6:8:end,:).'; C_all(:,6)=C6(:);
C7= C0 (7:8:end,:).'; C_all(:,7)=C7(:);
C8= C0 (8:8:end,:).'; C_all(:,8)=C8(:); 

C1_cmplex=C_all(1:2:end,:)+1i*C_all(2:2:end,:);


figure(1)
% plot(C_all(:,3));
plot(real(C1_cmplex(:,3)));hold on;
% plot(real(C1_cmplex(:,4)));

plot(C1(:,5));
plot(C1(:,6));
plot(real(C1_cmplex(:,7)));
plot(real(C1_cmplex(:,8)));
hold on;
xlabel('Sample index');
ylabel('ADC output')
xlim([1 200])


% min(A)/min(B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fre_sample=2.9184e9;
fre_sample=2.94912e9;
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

% figure(2)
% plot(fre/1e9,mag2db(abs(freA)));hold on;
% plot(fre/1e9,mag2db(abs(freB)));
% plot(fre/1e9,mag2db(abs(freC)));


xlabel('Frequency (GHz)');
grid on;

[Y,I]=max(abs(freD));
% I=270;
phase_1(frame_ind)=angle(freB(I)/freD(I))/pi*180; %EVB 2 vs 4 i.e phase difference of antenna 5 relative to antenna 6a
phase_2(frame_ind)=angle(freC(I)/freD(I))/pi*180; %EVB 3 vs 4 i.e phase difference of antenna 4 relative to antenna 6a
phase_3(frame_ind)=angle(freE(I)/freH(I))/pi*180; %EVB 5 vs 8 i.e phase difference of antenna 1 relative to antenna 6b
phase_4(frame_ind)=angle(freF(I)/freH(I))/pi*180; %EVB 6 vs 8 i.e phase difference of antenna 2 relative to antenna 6b
phase_5(frame_ind)=angle(freG(I)/freH(I))/pi*180; %EVB 7 vs 8 i.e phase difference of antenna 3 relative to antenna 6b
phase_6(frame_ind)=angle(freD(I)/freH(I))/pi*180; %EVB 4 vs 8 i.e phase difference of antenna 6a relative to antenna 6b



phase_7(frame_ind)=angle(freE(I)/freG(I))/pi*180;
phase_8(frame_ind)=angle(freE(I)/freH(I))/pi*180;

% if frame_ind==10
% figure(20)
%     hold on;grid on;
%     % plot(fre/1e9,mag2db(abs(freB)));
%     plot(fre/1e9,mag2db(abs(freC)));
%     plot(fre/1e9,mag2db(abs(freH)));
%     plot(fre/1e9,mag2db(abs(freG)));
%     % plot(fre/1e9,mag2db(abs(freF)));
%     % plot(fre/1e9,mag2db(abs(freG)));
%     % plot(fre/1e9,mag2db(abs(freH)));
%     xlabel('Frequeny (GHz)')
% end

end

%Y
P_diff1=(mod(phase_1+180,360)-180).';
P_diff2=(mod(phase_2+180,360)-180).';
P_diff3=(mod(phase_3+180,360)-180).';
P_diff4=(mod(phase_4+180,360)-180).';

P_diff5=(mod(phase_5+180,360)-180).';
P_diff6=(mod(phase_6+180,360)-180).';
P_diff7=(mod(phase_7+180,360)-180).';
P_diff8=(mod(phase_8+180,360)-180).';

phase_diff_all(:,ang_ind+1+FOV/2)=[median(P_diff1) median(P_diff2) median(P_diff3) median(P_diff4) median(P_diff5) median(P_diff6) median(P_diff7) median(P_diff8)];

figure(3)
plot([1:512]+(ang_ind+FOV/2)*512,P_diff1,'-o');hold on;grid on;
plot([1:512]+(ang_ind+FOV/2)*512,P_diff2,'-o');
plot([1:512]+(ang_ind+FOV/2)*512,P_diff4,'-o');
% plot(P_diff7,'-o');
% plot(P_diff8,'-o');
xlabel('Frame number')
ylabel('Phase difference (deg)')

end

figure(4)
plot(start_angle:end_angle,unwrap(phase_diff_all(1,FOV/2+1+(start_angle:end_angle))),'-o');hold on;grid on;
plot(start_angle:end_angle,unwrap(phase_diff_all(2,FOV/2+1+(start_angle:end_angle))),'-o');
plot(start_angle:end_angle,unwrap(phase_diff_all(3,FOV/2+1+(start_angle:end_angle))),'-o');
plot(start_angle:end_angle,unwrap(phase_diff_all(4,FOV/2+1+(start_angle:end_angle))),'-o');
plot(start_angle:end_angle,unwrap(phase_diff_all(5,FOV/2+1+(start_angle:end_angle))),'-o');
plot(start_angle:end_angle,unwrap(phase_diff_all(6,FOV/2+1+(start_angle:end_angle))),'-o');
xlabel('Angle (deg)')
ylabel('Phase difference (deg)')
legend('Phase diff 1','Phase diff 2','Phase diff 3','Phase diff 4','Phase diff 5','Between tile')

% t=[0:1023];
% C=cos(512*t/1024*2*pi);
% 
% figure(2)
% hold on;grid on;
% plot(fre/1e9,mag2db(abs(freB)));
% plot(fre/1e9,mag2db(abs(freC)));
% plot(fre/1e9,mag2db(abs(freD)));
% plot(fre/1e9,mag2db(abs(freG)));
% plot(fre/1e9,mag2db(abs(freF)));
% plot(fre/1e9,mag2db(abs(freG)));
% plot(fre/1e9,mag2db(abs(freH)));
% xlabel('Frequeny (GHz)')
% ylabel('Magnitude (dB)')

n = length(fre);
figure(2)
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
