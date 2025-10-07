clear all;close all;clc;

path = 'U:\Direction_Finding\20250930_MaranaCalibrationLibrary_915MHz_360AZ_66_to_-6EL\'; %keep trailing backslash
% path ='.\Data\';  % Active path (local relative folder)
bin_files = dir(fullfile(path, '*.BIN'));
num_files = numel(bin_files);
% ================= LOAD CALIBRATION LIBRARY =================
load('0.915GHz_cached_library_data.mat')
frequency_lib = 915;

% ================= ALGORITHM SELECTION =================
algorithm_select='MUSIC';
% algorithm_select='CAPON';
% algorithm_select='Correlation';

% ================= ROI FREQUENCY INDEX SELECTION =================
% Different candidate frequency bins of interest
% roi_freq=1:1024;
% roi_freq=741:751;%2.45 GHz in 20250626_MaranaTest_AZ360_EL66_Step3_withLens_withEVB_2.456GHz_TestData_skipfirsttwo
% roi_freq=320:344;%2.36 GHz LTE
% roi_freq=604:643; % most of drone data 
roi_freq = 570:580; %915 MHz

% roi_freq=619;
% roi_freq=680:730;
% roi_freq=327;

% ================= FLAGS =================
plot_flag=0;      % Enable plotting
print_detail=0;   % Print debug/info
save_flag=1;      % Save results

% ================= PARAMETERS =================
N_frame=1024;     % Number of frames
id_name=1;       % File index

tic
% ================= MAIN LOOP =================
for ang_ind=1:num_files  % Placeholder loop (expandable if scanning angles)
    % ---- Read binary input file ----
    id_name
    fileID = fopen([path sprintf('%04d', id_name) '.BIN'], 'r', 'ieee-le');


    DMA=2; % DMA = 1 for new firmware, DMA = 0 for old firmware
    if DMA == 1
        % Read raw IQ data

        C = fread(fileID, Inf, 'int16');fclose(fileID);
        C0 = reshape(C,[8,length(C)/8]).'; % Reshape into 8 channels
        L_C0=length(C0);
    
        % Split into 8 channels (DMA ordering)
        C1= C0 (1:2:L_C0/4,:).'; C_all(:,1)=C1(:);
        C2= C0 (2:2:L_C0/4,:).'; C_all(:,2)=C2(:);
        C3= C0 (L_C0/4+1:2:L_C0/4*2,:).'; C_all(:,3)=C3(:);
        C4= C0 (L_C0/4+2:2:L_C0/4*2,:).'; C_all(:,4)=C4(:);
        C5= C0 (L_C0/2+1:2:L_C0/4*3,:).'; C_all(:,5)=C5(:);
        C6= C0 (L_C0/2+2:2:L_C0/4*3,:).'; C_all(:,6)=C6(:);
        C7= C0 (L_C0/4*3+1:2:L_C0,:).'; C_all(:,7)=C7(:);
        C8= C0 (L_C0/4*3+2:2:L_C0,:).'; C_all(:,8)=C8(:); 
    
        % Reconstruct complex samples
        C1_cmplex=C_all(1:2:end,:)+1i*C_all(2:2:end,:);
    elseif DMA == 0
        % Alternative format
        C = fread(fileID, Inf, 'int16');fclose(fileID);
        C0 = reshape(C,[8,length(C)/8]).';
        
        % Assign channels sequentially
        C1= C0 (1:8:end,:).'; C_all(:,1)=C1(:);
        C2= C0 (2:8:end,:).'; C_all(:,2)=C2(:);
        C3= C0 (3:8:end,:).'; C_all(:,3)=C3(:);
        C4= C0 (4:8:end,:).'; C_all(:,4)=C4(:);
        C5= C0 (5:8:end,:).'; C_all(:,5)=C5(:);
        C6= C0 (6:8:end,:).'; C_all(:,6)=C6(:);
        C7= C0 (7:8:end,:).'; C_all(:,7)=C7(:);
        C8= C0 (8:8:end,:).'; C_all(:,8)=C8(:); 
        
        % Reconstruct complex samples
        C1_cmplex=C_all(1:2:end,:)+1i*C_all(2:2:end,:);
    elseif DMA == 2
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
    end
  

    % ================= FREQUENCY AXIS =================
    fre_sample=2.94912e9;        % Original sample rate
    fre_sample=2.94912e9/12;     % Decimated sample rate
    fre=[0:1023]/1024*fre_sample; % Frequency vector


    % ================= FFT PROCESSING =================
    FFT_full=zeros(N_frame,1024,6);      % FFT results [frame, freq, RX]
    threshold_FFT=zeros(N_frame,1024,6); % Thresholding results
   
    for frame_ind=1:N_frame
        % Extract FFT per channel (mapping order based on antenna channels)
        % freA=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),1));
        freB=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),2)); % channel 5
        freC=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),3)); % channel 4
        freD=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),4)); % channel 6
        freE=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),5)); % channel 1
        freF=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),6)); % channel 2
        freG=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),7)); % channel 3
        % freH=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),8));
        
        % Reorder into 6 RX channels
        FFT_full(frame_ind,end:-1:1,1)=fftshift(freE);
        FFT_full(frame_ind,end:-1:1,2)=fftshift(freF);
        FFT_full(frame_ind,end:-1:1,3)=fftshift(freG);
        FFT_full(frame_ind,end:-1:1,4)=fftshift(freC);
        FFT_full(frame_ind,end:-1:1,5)=fftshift(freB);
        FFT_full(frame_ind,end:-1:1,6)=fftshift(freD);      
    end


    % ================= PLOT SPECTROGRAM =================
    
    Nrx=6;   % Number of RX channels used
    time_axis=[0:1023]*1/fre_sample*1e6;  % Time axis (us)
    fre_axis=fre/1e9+.900 - 245.76/2e3;               % Frequency axis (GHz), with LO offset
    

if(plot_flag)    
    for rx=1:Nrx
        test_STFT=abs(FFT_full(:,:,rx));  % Magnitude spectrogram for RX
        
        
            figure(98);
            subplot(2, 3, rx);
            imagesc(time_axis,fre_axis,mag2db(test_STFT.')); % Plot spectrogram
            set(gca, 'YDir', 'normal');  % Flip Y axis
            grid on;
            
            xlabel('Time (us)');
            ylabel('Frequency (GHz)');
            cb = colorbar;               
            ylabel(cb, 'Magnitude (dB)');    
            caxis([20 90]);              % Dynamic range for plotting
            title(sprintf('Spectrogram of RX %d', rx));
    end
end
     
    % ================= CFAR DETECTION =================
    % n = length(fre);
    % frame_peak_indices_all = cell(1, N_frame);
    % frame_all_indices=[];
 
    offset=10;

    mag=abs(FFT_full);
    noise_all0=median(mag,2);
    noise_all=repmat(noise_all0,[1 1024 1]);

    threshold_cfar = offset * noise_all;
    threshold_FFT = mag > threshold_cfar;


    % 
    % for frame_idx=1:N_frame
    %     peak_indices_all = cell(1, Nrx);
    % 
    %     for rx = 1:Nrx
    % 
    %         test=FFT_full(frame_idx,:,rx);   % FFT slice per frame & RX
    %         mag=abs(test);                  % Magnitude spectrum
    % 
    %         CA_CFAR_detection_median();            % Call CFAR function (user-defined)
    % 
    %         detected_indices = find(cfar_output == 1); % Get detected peaks
    %         peak_indices_all{rx} = detected_indices;                        
    %     end
    % 
    %     % Mark detected peaks in threshold matrix
    %     for rx = 1:Nrx
    %         threshold_FFT(frame_idx,peak_indices_all{rx},rx)=1;
    %     end
    % end

    % ================= COMBINE DETECTIONS ACROSS RX =================
    final_threshold=ones(N_frame,1024);
    for rx=1:Nrx
        test_threshold=threshold_FFT(:,:,rx);
        
        if(plot_flag)
            figure(97);
            subplot(2, 3, rx);
            imagesc(time_axis,fre_axis,(test_threshold.')); % Binary detections
            set(gca, 'YDir', 'normal');
            grid on;
            xlabel('Time (us)');
            ylabel('Frequency (GHz)');
            colormap(flipud(gray))
            cb = colorbar;              
            title(sprintf('Binary threshold of RX %d', rx));
        end
        
        final_threshold=final_threshold.*test_threshold; % Common detections
    end
    
    % Final binary threshold (intersection across RXs)
    if(plot_flag)
        figure()
        imagesc(time_axis,fre_axis,(final_threshold.'));
        set(gca, 'YDir', 'normal');
        grid on;
        xlabel('Time (us)');
        ylabel('Frequency (GHz)');
        colormap(flipud(gray))
        cb = colorbar;              
        title(sprintf('Final threshold'));
    end


    % ================= ANGLE SCAN GRID =================
    if(frequency_lib == 2456)
        AZ_start = 180; AZ_end = -180; AZ_step = -3;
        EL_start = 66; EL_end = 0; EL_step = -3;
    elseif(frequency_lib ==2360 | frequency_lib == 915)
        AZ_start = 180; AZ_end = -180; AZ_step = -3;
        EL_start = 66; EL_end = -6; EL_step = -3;
    end

    % Build azimuth and elevation scan tables
    AZ_data = AZ_start:AZ_step:AZ_end;
    AZ_steps = length(AZ_data);
    EL_data = EL_start:EL_step:EL_end;
    EL_steps = length(EL_data);
    
    AZ_table=AZ_data;
    EL_table=EL_data;
    length_el=length(EL_table);
    length_az=length(AZ_table);
    
    

    
    % ================= INITIALIZE RESULT ARRAYS =================
    AZ_result=[];
    AZ_result_itp=[];
    EL_result=[];
    EL_result_itp=[];
    summary_result=[];

    % ================= MAIN PEAK PROCESSING LOOP =================
    for peak_idx = 1:length(roi_freq)
        peak = roi_freq(peak_idx);
        fre_peak=fre(peak)/1e9+0.915 - 245.76/2e3; % Convert bin index → GHz

        % Find frames where CFAR detected the peak
        frames = find(final_threshold(:,peak) ==1);
        
        if(print_detail)
            fprintf('Peak %d is found in frame(s): %s\n', peak, mat2str(frames'));
        end
        
        if (length(frames)>1)
            % Extract complex FFT data for selected frames
            test_data=squeeze(FFT_full(frames,peak,:));
            X=test_data.'; 
            R = (X * X') / length(frames); % Covariance matrix
            
            % ================= DOA ALGORITHM =================
            if strcmp(algorithm_select,'MUSIC')
                % --------- MUSIC Algorithm ---------
                [Evecs, Evals] = eig(R); % Eigen decomposition
                [Evals_sorted, idx] = sort(diag(Evals), 'descend');
                Evecs_sorted = Evecs(:, idx);
                d=2; % number of signals assumed
                En = Evecs_sorted(:, d+1:end); % Noise subspace
                
                % MUSIC spectrum evaluation
                Pmusic = zeros(length_el, length_az);
                for el_idx = 1:length_el
                    for az_idx = 1:length_az
                        steering = Lib_Complex(:, az_idx, el_idx);  
                        Pmusic(el_idx, az_idx) = 1 / abs(steering' * (En * En') * steering);
                    end
                end
                
                Pmusic_dB = mag2db(abs(Pmusic));
                abs_spectrum=abs(Pmusic);
    
                % Plot MUSIC spectrum if only 1 ROI frequency
                if(length(roi_freq)==1)
                    figure;
                    imagesc(AZ_table, EL_table, Pmusic_dB);
                    set(gca, 'YDir', 'normal');
                    xlabel('Azimuth (°)'); ylabel('Elevation (°)');
                    cb = colorbar; ylabel(cb, 'Spatial Spectrum (dB)');
                    grid on; title('2D MUSIC Spectrum');
                    
                    [AZ_grid, EL_grid] = meshgrid(AZ_table, EL_table);
                    figure; surf(AZ_grid, EL_grid, Pmusic_dB, 'EdgeColor', 'none');
                    xlabel('Azimuth (°)'); ylabel('Elevation (°)'); zlabel('Spatial Spectrum (dB)');
                    title('2D MUSIC Spectrum (3D View)');
                    cb = colorbar; ylabel(cb, 'Spatial Spectrum (dB)');
                    view(45, 45); shading interp;
                end
                
            elseif strcmp(algorithm_select,'CAPON')
                % --------- CAPON Algorithm ---------
                try
                    R_inv = pinv(R);  % Pseudo-inverse covariance matrix
                catch
                    fprintf('ERROR: Singular matrix\n');
                    ADSINR = -3;  
                    return;
                end
               
                PAD = zeros(length_el, length_az);
                for el_idx = 1:length_el
                    for az_idx = 1:length_az
                        steering = Lib_Complex(:, az_idx, el_idx);
                        PAD(el_idx, az_idx) = (steering' * R_inv) * steering; 
                    end
                end
                PCapon=1./PAD;
                
                abs_spectrum=abs(PCapon);
                PCapon_dB = mag2db(abs(PCapon));
                
                if(length(roi_freq)==1)
                    figure;
                    imagesc(AZ_table, EL_table, PCapon_dB);
                    set(gca, 'YDir', 'normal');
                    xlabel('Azimuth (°)'); ylabel('Elevation (°)');
                    cb = colorbar; ylabel(cb, 'Spatial Spectrum (dB)');
                    grid on; title('2D CAPON Spectrum');
                    
                    [AZ_grid, EL_grid] = meshgrid(AZ_table, EL_table);
                    figure; surf(AZ_grid, EL_grid, PCapon_dB, 'EdgeColor', 'none');
                    xlabel('Azimuth (°)'); ylabel('Elevation (°)'); zlabel('Spatial Spectrum (dB)');
                    title('2D CAPON Spectrum (3D View)');
                    cb = colorbar; ylabel(cb, 'Spatial Spectrum (dB)');
                    view(45, 45); shading interp;
                end
                
            elseif strcmp(algorithm_select,'Correlation')
                % --------- Correlation-based DOA ---------
                rxSignal2D_FFT_rx1=test_data(:,1);
                rxSignal2D_FFT_test=test_data./rxSignal2D_FFT_rx1;
                rxSignal2D_FFT_test_mean=mean(rxSignal2D_FFT_test,1).';
                weighting=[1 1 1 1 1 1]';
                
                for az=1:AZ_steps
                    for el=1:EL_steps
                        r=corrcoef(rxSignal2D_FFT_test_mean.*weighting, ...
                                   squeeze(Lib_Complex(:,az,el)).*weighting);
                        coe(az,el)=r(1,2);
                    end
                end
                
                abs_spectrum=reshape(abs(coe),[length_az length_el]);
                
                if(length(roi_freq)==1)
                    figure();
                    imagesc(EL_table, AZ_table, mag2db(abs_spectrum));
                    set(gca, 'YDir', 'normal');
                    xlabel('Elevation (deg)'); ylabel('Azimuth (deg)');
                    colorbar; colormap jet;
                    title('Color Field Map');
                    grid on;
                end
                abs_spectrum=abs_spectrum.';
            end
           
            % ================= PEAK SEARCH & INTERPOLATION =================
            [max_val, max_idx] = max(abs_spectrum(:));
            [peak_el, peak_az] = ind2sub(size(abs_spectrum), max_idx);
            
            AZ_peak = AZ_table(peak_az);
            EL_peak = EL_table(peak_el);
            
            cf_reshape = abs_spectrum;       
            AZ_step = AZ_table(2) - AZ_table(1);
            EL_step = EL_table(2) - EL_table(1);
            
            % Interpolation in Azimuth
            delta_az = 0;
            if peak_az > 1 && peak_az < length_az
                a = cf_reshape(peak_el, peak_az - 1);
                b = cf_reshape(peak_el, peak_az);
                c = cf_reshape(peak_el, peak_az + 1);
                if a ~= c
                    delta_az = 0.5 * (a - c) / (a - 2*b + c);
                end
            end
            
            % Interpolation in Elevation
            delta_el = 0;
            if peak_el > 1 && peak_el < length_el
                a1 = cf_reshape(peak_el - 1, peak_az);
                b1 = cf_reshape(peak_el, peak_az);
                c1 = cf_reshape(peak_el + 1, peak_az);
                if a1 ~= c1
                    delta_el = 0.5 * (a1 - c1) / (a1 - 2*b1 + c1);
                end
            end
            
            % Final interpolated peaks
            AZ_peak_itp = AZ_table(peak_az) + delta_az * AZ_step;
            EL_peak_itp = EL_table(peak_el) + delta_el * EL_step;
                        
            % Store results
            AZ_result=[AZ_result;AZ_peak];
            AZ_result_itp=[AZ_result_itp;AZ_peak_itp];
            EL_result=[EL_result;EL_peak];
            EL_result_itp=[EL_result_itp;EL_peak_itp];
        
            summary_result=[summary_result;peak fre_peak length(frames) AZ_peak EL_peak AZ_peak_itp EL_peak_itp];

        end
    end

    toc

    if ~isempty(AZ_result)
        final_AZ=median(AZ_result);
        final_AZ_ITP=median(AZ_result_itp)
    
        final_EL=median(EL_result);
        final_EL_ITP=median(EL_result_itp);
        if(print_detail)
            fprintf('Final Results:\n');
            fprintf('  AZ      = %.3f\n', final_AZ);
            fprintf('  AZ ITP  = %.3f\n', final_AZ_ITP);
            fprintf('  EL      = %.3f\n', final_EL);
            fprintf('  EL ITP  = %.3f\n', final_EL_ITP); 
        end
    end
    if(save_flag==1)
        if(frequency_lib==2456)
            save_folder='Result_2p45_lib';
            save_folder = [save_folder '_' algorithm_select];
            if ~exist(save_folder, 'dir')
                mkdir(save_folder);
            end
            save([save_folder '/' num2str(id_name) '_2p45_result.mat'], 'summary_result','final_AZ','final_AZ_ITP','final_EL','final_EL_ITP')
            id_name=id_name+1;
        elseif(frequency_lib==2360)
            save_folder='Result_2p36_lib';
            save_folder = [save_folder '_' algorithm_select];
            if ~exist(save_folder, 'dir')
                mkdir(save_folder);
            end
            save([save_folder '/' num2str(id_name) '_2p36_result.mat'], 'summary_result','final_AZ','final_AZ_ITP','final_EL','final_EL_ITP')
            id_name=id_name+1;
        elseif(frequency_lib == 915)
            save_folder='Result_0p915_lib';
            save_folder = [save_folder '_' algorithm_select];
            if ~exist(save_folder, 'dir')
                mkdir(save_folder);
            end
            save([save_folder '/' num2str(id_name) '_0p915_result.mat'], 'summary_result','final_AZ','final_AZ_ITP','final_EL','final_EL_ITP')
            id_name=id_name+1;
        end
    end
    
    if ~isempty(AZ_result)
    
        distance=1;
        X_plot=distance*cosd(final_AZ_ITP);Y_plot=distance*sind(final_AZ_ITP);
        X_plot2=distance*cosd(final_EL_ITP);Y_plot2=distance*sind(final_EL_ITP);
    
        figure(1)
        set(gcf, 'Position',  [100, 100, 1200, 500]);
    
        subplot(1,2,1)
        hold on;
        plot(X_plot,Y_plot,'xk','MarkerSize',12,'Linewidth',2);
        plot([0,X_plot],[0,Y_plot],'r:','Linewidth',2);hold off;
        grid on;
        xlim([-1.5 1.5]);
        ylim([-1.5 1.5]);
        xlabel('X (m)');ylabel('Y (m)');
        subplot(1,2,2)
        hold on;
        plot(X_plot2,Y_plot2,'xk','MarkerSize',12,'Linewidth',2);
        plot([0,X_plot2],[0,Y_plot2],'r:','Linewidth',2);hold off;
        grid on;
        xlim([-1.5 1.5]);
        ylim([-1.5 1.5]);
        xlabel('X (m)');ylabel('Z (m)');
    
        fontsize(12, "points");
    end

end



