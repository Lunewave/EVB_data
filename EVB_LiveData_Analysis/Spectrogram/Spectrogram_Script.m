%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all;close all;clc;
path = 'U:\Falcon_Project\20250714_LWOfficeTest_AZ20_EL0_Step5_withLens_withEVB_DroneDebug_wDMA\Room_with_drone_and_remote\';
frame = 4;

DMA=1;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen([path num2str(frame,'%04d') '.BIN'], 'r', 'ieee-le');
if DMA
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
else

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
end

antenna_order = [5 4 6 1 2 3];


% Parameters for spectrogram
window_size = 1024;         % Length of FFT window
overlap = round(0.75 * window_size);  % 75% overlap
nfft = 1024;                % Number of FFT points
fs = 245.76e6;  % 245.76 MHz

% Loop over each channel
fc = 2.4e9;  % example center frequency

for ch = 2:7
    subplot(2, 3, ch-1);

    signal = C1_cmplex(:, ch);

    [s, f, t] = spectrogram(signal, window_size, overlap, nfft, fs, 'centered');
    rf_freq = f + fc;

    % Flip vertically so low freq is at bottom
    imagesc(t, rf_freq(end:-1:1)/1e9, 20*log10(abs(s(end:-1:1, :))));
    axis xy;
    title(['Spectrogram - Antenna ' num2str(antenna_order(ch-1))]);
    xlabel('Time (s)');
    ylabel('Frequency (GHz)');
    colormap turbo;
    colorbar;
end


sgtitle('Spectrograms of Antennas 1–6 (RF-Aligned)');

