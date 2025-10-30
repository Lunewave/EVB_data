%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all;close all;clc;
path = 'U:\Direction_Finding\20250930_MaranaCalibrationLibrary_915MHz_360AZ_66_to_-6EL\';
frame = 2000;

num_frames = 0;


for i = frame:frame+num_frames

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fileID = fopen([path num2str(i,'%04d') '.BIN'], 'r', 'ieee-le');

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

    antenna_order = [5 4 6 1 2 3];
    figure(i+1000)
        % Parameters for spectrogram
        window_size = 1024;         % Length of FFT window
        overlap = round(0.75 * window_size);  % 75% overlap
        nfft = 1024;                % Number of FFT points
        fs = 245.76e6;  % 245.76 MHz
            
        fc = 0.9e9;  % or whatever center frequency you're using
        
    for ch = 2:7
        subplot(2, 3, ch-1);
    
        signal = conj(C1_cmplex(:, ch));
    
        [s, f, ~] = spectrogram(signal, window_size, overlap, nfft, fs, 'centered');
        rf_freq = f + fc;
    
        % Construct custom x-axis: one label per column in spectrogram
        num_frames = size(s, 2);
        frame_numbers = linspace(1, 1024, num_frames);  % span full 1024 frames
    
        % Plot spectrogram with RF-aligned y-axis and frame-based x-axis
        imagesc(frame_numbers, rf_freq(end:-1:1)/1e9, 20*log10(abs(s(end:-1:1, :))));
        axis xy;
        title(['Spectrogram - Antenna ' num2str(antenna_order(ch-1))]);
        xlabel('Frame Number');
        ylabel('Frequency (GHz)');
        colorbar;
        caxis([20 90]);              % Dynamic range for plotting
    end
    
    sgtitle('Spectrograms of Antennas 1–6 (RF-Aligned by Frame)');
    fre=[0:1023]/1024*fs;
    a = fre/1e9+(fc - fs/2)/1e9;
    n = length(fre);
    a = a(end:-1:1);
    a = [a(n/2+1:end),a(1:n/2)];
    [a_sorted, idx] = sort(a, 'ascend');


    for frame_ind=1:1024
        % freA=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),1));
        freB=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),2)); %CH2
        freC=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),3)); %CH3
        freD=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),4)); %CH4
        freE=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),5)); %CH5
        freF=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),6)); %CH6
        freG=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),7)); %CH7
        % freH=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),8)); %CH8
    end
    figure(2000+i)
    hold on; grid on;
    plot(a_sorted, mag2db(abs(freB(idx))));
    plot(a_sorted, mag2db(abs(freC(idx))));
    plot(a_sorted, mag2db(abs(freD(idx))));
    plot(a_sorted, mag2db(abs(freE(idx))));
    plot(a_sorted, mag2db(abs(freF(idx))));
    plot(a_sorted, mag2db(abs(freG(idx))));
    % yline(65, '--r', 'Threshold');
    xlabel('Frequency (GHz)')
    ylabel('Magnitude (dB)')
    title(['FFT of file: ' num2str(i)]);
    legend('Antenna 5','Antenna 4', 'Antenna 6', 'Antenna 1', 'Antenna 2', 'Antenna 3', 'Location', 'best');
end