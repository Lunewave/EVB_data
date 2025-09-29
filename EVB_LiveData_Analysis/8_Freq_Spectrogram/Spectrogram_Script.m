%%%%%%%%%%%%%%%%%CHANGE THESE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all;close all;clc;
path = 'U:\Direction_Finding\20250926_udp_output_ConferenceRoomCollection\';
frame = 1; %Starting frame
num_frames_to_plot = 8; %how many frames to make spectrograms for







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bin_files = dir(fullfile(path, '*.BIN'));
file1 = bin_files(1).name;
file1 = strsplit(file1, '.');
file1 = str2num(file1{1});


fre_list = [2.4 2.3 2.2 2 1.8 1.6 1.3 0.9]* 1e9;
fs = 2.94912e9;
bw = fs/12; %245.76 MHz
nyquist = fs/2;



for i = frame:frame+num_frames_to_plot

   
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
    fre_idx = mod(i - file1, 8) + 1;
    fc = fre_list(fre_idx);
    for ch = 2:7
        subplot(2, 3, ch-1);
        if fc<nyquist
            signal = conj(C1_cmplex(:, ch));
        else
            signal = C1_cmplex(:, ch);
        end
    
        [s, f, ~] = spectrogram(signal, window_size, overlap, nfft, bw, 'centered');
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
        colormap turbo;
        colorbar;
    end
    
    sgtitle(['Spectrograms of Antennas 1–6 (RF-Aligned by Frame) ' num2str(fre_list(fre_idx)/1e9) ' GHz']);

    
 
end