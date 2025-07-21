%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all;close all;clc;
path = 'U:\Falcon_Project\20250711_MaranaTest_AZ360_EL0_Step5_withLens_withEVB_2.456GHz_DroneTest_r-10_h-10\';
frame = 37;

DMA=1;
num_frames = 0;


for i = frame:frame+num_frames

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fileID = fopen([path num2str(i,'%04d') '.BIN'], 'r', 'ieee-le');
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
    figure(i+1000)
    if DMA
        % Parameters for spectrogram
        window_size = 1024;         % Length of FFT window
        overlap = round(0.75 * window_size);  % 75% overlap
        nfft = 1024;                % Number of FFT points
        fs = 245.76e6;  % 245.76 MHz
            
        fc = 2.4e9;  % or whatever center frequency you're using
        
        for ch = 2:7
            subplot(2, 3, ch-1);
        
            signal = C1_cmplex(:, ch);
        
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
            colormap turbo;
            colorbar;
        end
        
        sgtitle('Spectrograms of Antennas 1–6 (RF-Aligned by Frame)');

    
    else
        fc = 2.4e9;                     % Center frequency (adjust as needed)
        fs = 245.76e6;                  % Sampling rate
        window_size = 1024;
        overlap = 0;                   % No overlap to avoid smearing across bursts
        nfft = 1024;
        samples_per_frame = 1024;
        frames_per_burst = 8;
        samples_per_burst = samples_per_frame * frames_per_burst;
        num_bursts = 1024 / frames_per_burst;
        antenna_order = [5 4 6 1 2 3];
        
        for ch = 2:7
            subplot(2, 3, ch - 1);
            full_s = [];
            full_x = [];
        
            for burst = 1:num_bursts
                idx_start = (burst - 1) * samples_per_burst + 1;
                idx_end   = burst * samples_per_burst;
                burst_signal = C1_cmplex(idx_start:idx_end, ch);
        
                [s, f, ~] = spectrogram(burst_signal, window_size, overlap, nfft, fs, 'centered');
        
                % Append results
                full_s = [full_s, s];
        
                % Map time columns of this burst to frame numbers
                burst_frames = ((burst - 1) * frames_per_burst) + linspace(1, frames_per_burst, size(s, 2));
                full_x = [full_x, burst_frames];
            end
        
            rf_freq = f + fc;
        
            % Plot the combined spectrogram
            imagesc(full_x, rf_freq(end:-1:1)/1e9, ...
                    20*log10(abs(full_s(end:-1:1, :))));
            axis xy;
            title(['Spectrogram - Antenna ' num2str(antenna_order(ch - 1))]);
            xlabel('Frame Number');
            ylabel('Frequency (GHz)');
            % Force Y-axis ticks every 0.01 GHz
            ytick_vals = round(min(rf_freq)/1e6)*1e6 : 10e6 : round(max(rf_freq)/1e6)*1e6;
            yticks(ytick_vals / 1e9);
            colormap turbo;
            colorbar;
        
            hold on;
            for k = 8:8:1024
                xline(k + 0.5, 'k-', 'LineWidth', 0.5);
            end
        end
        
        sgtitle('Burst-wise Spectrograms of Antennas 1–6 (No Overlap)');
    end
end