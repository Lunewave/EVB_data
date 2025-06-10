clear all;close all;clc;

fileID = fopen('Old data\ADC_DATA_05052025_1.BIN', 'r', 'ieee-le');

A = fread(fileID, Inf, 'int16');fclose(fileID);

% figure(1)
% plot(A);
% hold on;

% for file_ind=0:7
%     fileID = fopen(['ADC_' num2str(file_ind)  '.BIN'], 'r', 'ieee-le');
%     % fileID = fopen(['ADC_DATA_05052025_1.BIN'], 'r', 'ieee-le');
% 
%     B = fread(fileID, Inf, 'int16');fclose(fileID);
% 
%     figure(1)
%     plot(B);
% 
% end
% grid on;

% fileID = fopen('ADC.BIN', 'r', 'ieee-le');
% fileID = fopen('DATA_R6\ADC.BIN', 'r', 'ieee-le');

for ang_ind=0:10
fileID = fopen(['DATA_49_7ch_drone_2p4G_image2\' num2str(ang_ind,'%04d') '.BIN'], 'r', 'ieee-le');
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

% plot(C1(:,5));
% plot(C1(:,6));
% plot(real(C1_cmplex(:,7)));
% plot(real(C1_cmplex(:,8)));
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
I=270;
phase_1(frame_ind)=angle(freB(I)/freD(I))/pi*180;
phase_2(frame_ind)=angle(freC(I)/freD(I))/pi*180;
phase_3(frame_ind)=angle(freE(I)/freH(I))/pi*180;
phase_4(frame_ind)=angle(freF(I)/freH(I))/pi*180;
% phase_4(frame_ind)=angle(freC(I))/pi*180;

phase_5(frame_ind)=angle(freG(I)/freH(I))/pi*180;
phase_6(frame_ind)=angle(freD(I)/freH(I))/pi*180;
phase_7(frame_ind)=angle(freE(I)/freG(I))/pi*180;
phase_8(frame_ind)=angle(freE(I)/freH(I))/pi*180;

if frame_ind==10
figure(20)
    hold on;grid on;
    % plot(fre/1e9,mag2db(abs(freB)));
    plot(fre/1e9,mag2db(abs(freC)));
    plot(fre/1e9,mag2db(abs(freH)));
    plot(fre/1e9,mag2db(abs(freG)));
    % plot(fre/1e9,mag2db(abs(freF)));
    % plot(fre/1e9,mag2db(abs(freG)));
    % plot(fre/1e9,mag2db(abs(freH)));
    xlabel('Frequeny (GHz)')
end

end

%Y
P_diff1=(mode(phase_1+180,360)-180).';
P_diff2=(mode(phase_2+180,360)-180).';
P_diff3=(mode(phase_3+180,360)-180).';
P_diff4=(mode(phase_4+180,360)-180).';

P_diff5=(mode(phase_5+180,360)-180).';
P_diff6=(mode(phase_6+180,360)-180).';
P_diff7=(mode(phase_7+180,360)-180).';
P_diff8=(mode(phase_8+180,360)-180).';

phase_diff_all(:,ang_ind+1)=[median(P_diff1) median(P_diff2) median(P_diff3) median(P_diff4)...
    median(P_diff5) median(P_diff6) median(P_diff7) median(P_diff8)];

figure(3)
plot([1:512]+ang_ind*512,P_diff1,'-o');hold on;grid on;
plot([1:512]+ang_ind*512,P_diff2,'-o');
plot([1:512]+ang_ind*512,P_diff4,'-o');
% plot(P_diff7,'-o');
% plot(P_diff8,'-o');
xlabel('Frame number')
ylabel('Phase difference (deg)')

end

figure(4)
plot(phase_diff_all(1,:),'-o');hold on;grid on;
plot(phase_diff_all(2,:),'-o');
plot(phase_diff_all(3,:),'-o');
plot(phase_diff_all(4,:),'-o');
plot(phase_diff_all(5,:),'-o');
plot(phase_diff_all(6,:),'-o');
xlabel('Ange (deg)')
ylabel('Phase difference (deg)')
legend('Phase diff 1','Phase diff 2','Phase diff 3','Phase diff 4','Phase diff 5','Between tile')

% t=[0:1023];
% C=cos(512*t/1024*2*pi);
% 
figure(2)
hold on;grid on;
% plot(fre/1e9,mag2db(abs(freB)));
plot(fre/1e9,mag2db(abs(freC)));
plot(fre/1e9,mag2db(abs(freD)));
plot(fre/1e9,mag2db(abs(freG)));
% plot(fre/1e9,mag2db(abs(freF)));
% plot(fre/1e9,mag2db(abs(freG)));
% plot(fre/1e9,mag2db(abs(freH)));
xlabel('Frequeny (GHz)')