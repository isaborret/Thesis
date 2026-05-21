clear all;close all;clc;

%% TIME SERIES EXTRACTION - INTRABEAT - RESTING STATE (CHIETI)
addpath(genpath("C:\Users\Isa\OneDrive\Documenten\Unif 2e master\Masterproef_II\code\preprocessing\FUNCTIONS"));

%% Parameters
% frequency bands
name_bands = {'delta', 'theta', 'alpha', 'sigma', 'beta'};
ybands = [0.5, 4.0, 8.0, 12.0, 16.0, 25.0]; 


% Computes the bandwidth of each band (e.g. delta = 4-0.5 = 3.5 Hz)
% This is used later to scale the power estimates correctly — multiplying mean PSD by bandwidth gives total power in that band.
for ibands = 1:length(ybands)-1 
    base_banda(ibands) = ybands(ibands+1)-ybands(ibands);
end

%base_tot =  ybands(end) - ybands(1); % Total bandwidth across all bands combined, used for normalisation

% Parameters for the weighted covariance PSD estimation
nfft=512; %512 frequency bins
Bw= 0.02; % bandwidth 0.02 Hz
window= 'hamming'; 

% HRV parameters
pfilter = 0.94;
spectrum = [0.05, 0.15, 0.4];
base_LF = spectrum(2)-spectrum(1);
base_HF = spectrum(3)-spectrum(2); % HF band width = 0.25 Hz
base_tot= spectrum(3)-spectrum(1); % total HRV band width

% Beat selection => check!!!!!!!! (new)
max_beats = 300; % conservative limit for 5-minute epochs


% sampling rate and channel labels taken from first subject (assumed identical across subjects)
% fs = data(1).fs;
% ch_label = data(1).ch_label;
% resampling = 0; % 0= resampling is off
% passo = 1; %% 1 s - 1 Hz, 2 s - 0.5 Hz, 0.5 s - 2 Hz

%% Load and process subjects
input_dir = 'C:\Users\Isa\OneDrive\Documenten\Unif 2e master\Masterproef_II\data\derivatives\preprocessed_eeg\clean_mat'
%input_dir = 'C:\Users\Isa\OneDrive\Documenten\Unif 2e master\Masterproef_II\data\derivatives\5_minuteExtractionRawEEGandECG';
files     = dir(fullfile(input_dir, 'sub*.mat'));
n_subj    = length(files);
% Subjects to exclude from processing
ignore_subs = {'sub001'};

%% loop
for isubj = 1:n_subj
    % extract subject number from filename
    subnum = files(isubj).name(4:6);
    
    % Skip excluded subjects
    if ismember(['sub' subnum], ignore_subs)
        fprintf('Skipping subject sub%s (%s) - excluded.\n', subnum, files(isubj).name);
        continue
    end

    fprintf('Processing subject %d of %d: %s\n', isubj, n_subj, files(isubj).name);
    

    %% Load subject data
    tmp_rr = load(fullfile(input_dir, files(isubj).name));

    t           = tmp_rr.t;
    ECG         = tmp_rr.ECG_data(:);    
    EEG         = tmp_rr.EEG_data;       
    EEG_labels  = tmp_rr.EEG_labels;
    fs          = tmp_rr.fs;             % 256 Hz

    % Find Cz channel index
    cz_idx = find(strcmpi(EEG_labels, 'Cz'));
    if isempty(cz_idx)
        warning('Cz not found for subject %s, skipping.', files(isubj).name);
        continue
    end

    % Find Pz channel index
    pz_idx = find(strcmpi(EEG_labels, 'Pz'));
    if isempty(pz_idx)
        warning('Pz not found for subject %s, skipping.', files(isubj).name);
        continue
    end


    %% Store metadata
    results(isubj).subject_id  = tmp_rr.subject_id;
    results(isubj).dream       = tmp_rr.dream;
    results(isubj).sleep_stage = tmp_rr.sleep_stage;
    
    %% ECG pipeline
    %%% filtering: 4th order Butterworth filter (passband: 0.1-20Hz) --> in questo dataset già filtrato
    % t = data(isubj).t;
    % ECG = data(isubj).ECG;

    % Step 1: R-peak detection (Pan-Tompkins)
    [t_PT, ECG_PT, posMassimoPT, tMassimoPT, massimoPT] = functECG_PanTompkins(t, ECG); % Runs Pan-Tompkins on raw ECG → candidate R-peak positions. stayed the same

    % Step 2: R-peak refinement
    [posRPeakECG, tPosRPeakECG, massimoRPeakECG] = functECG_findRPeak(t, ECG, posMassimoPT);

    % Step 3: RR interval extraction
    [beatsRRI,tSeriesRRI,seriesRRI]=funct_ExtractTimeSeries(tPosRPeakECG,0);
    mean_RR = mean(seriesRRI); % Computes mean RR for the ±25% threshold

    % Step 4: Valid beat selection (+-25% mean RR threshold)
    % Iterates through all detected beats (u), accepting each beat into sel only if its RR interval is within ±25% of the mean. 
    % Keeps going until 300 valid beats are collected (k). seriesRRI_sel and tSeriesRRI_sel accumulate the valid RR intervals and their timestamps.
    k = 1;
    u = 1;
    sel = [];

    while k <= max_beats && u <= length(seriesRRI)
        rr_val = seriesRRI(u);
        if (rr_val > 4/5*mean_RR) && (rr_val < 6/5*mean_RR) % changed it to the 20% threshold
            k = k+1;
            sel = [sel u];
            seriesRRI_sel = seriesRRI(sel);
            tSeriesRRI_sel = tSeriesRRI(sel);
        end
        u = u+1;  
    end

    % Check sufficient beats were found
    if length(sel) < 200
        warning('Subject %s: only %d valid beats found, skipping.', ...
            files(isubj).name, length(sel));
        continue
    end

    n_beats = length(sel);

    fprintf('Before removeOutlier: seriesRRI_sel=%d, tSeriesRRI_sel=%d, n_beats=%d\n', ...
    length(seriesRRI_sel), length(tSeriesRRI_sel), n_beats);

    % Step 5: Secondary outlier correction on RR series => new!!
    seriesRRI_sel_corrected = removeOutlier(seriesRRI_sel, [], 'linear');
    %check whether series is shorter than expected
    if length(seriesRRI_sel_corrected) < n_beats
        seriesRRI_sel_corrected(end+1:n_beats) = seriesRRI_sel_corrected(end); % creates indices from the current last position plus one up to n_beats
                                                                               % extends the series by one sample by repeating the last value
    end
    seriesRRI_sel = seriesRRI_sel_corrected; % Overwrites seriesRRI_sel with the padded version 

    % Store RR series
    results(isubj).RRI      = seriesRRI_sel;
    results(isubj).tRRI     = tSeriesRRI_sel;
    results(isubj).n_beats  = n_beats;

    %% HF-HRV extraction via SPWVD 
    % Step 6: AR high-pass filter RR series
    RRI_filt = AR_filter(seriesRRI_sel', 1, pfilter);  % high-pass output, removes slow drift

    % Effective sampling rate of RR series
    fs_rri = 1000 / mean(seriesRRI_sel);

    % Step 7: Smoothed pseudo-Wigner-Ville distribution
    [TFR_RR, f_rr] = wvd(hilbert(RRI_filt - mean(RRI_filt)), fs_rri, ...
        'smoothedPseudo', hamming(51), 'MinThreshold', 50);

    % Step 8: Integrate over HF band and normalise
    ind_HF  = find(f_rr >= 0.15 & f_rr < 0.4);
    ind_tot = find(f_rr >= 0.04 & f_rr < 0.4);

    P_RRHF = zeros(1, size(TFR_RR, 2));
    for l = 1:size(TFR_RR, 2)
        P_RRHF(l) = (mean(TFR_RR(ind_HF, l))^2) * base_HF ./ ...
                    (mean(TFR_RR(ind_tot, l))^2) * f_rr(end);
    end

    %results(isubj).HFHRV = P_RRHF;
    %% downsampling
    hfhrv = downsample(P_RRHF', 2)';
    results(isubj).HFHRV = hfhrv(1:n_beats);
    if length(results(isubj).HFHRV) ~= n_beats
        fprintf('WARNING: Subject %s HFHRV length %d does not match n_beats %d\n', files(isubj).name, length(results(isubj).HFHRV), n_beats);
        hfhrv = results(isubj).HFHRV; % Copies HFHRV into a temporary variable so it can be modified
        hfhrv(end+1:n_beats) = hfhrv(end); % extends the series to n_beats by repeating the last value
        results(isubj).HFHRV = hfhrv; % Writes the padded series back into the results struct
    end

    %% EEG band power extraction - Cz and Pz
    for ich = 1:2
        if ich == 1
            chan_idx  = cz_idx;
            chan_name = 'Cz';
        else
            chan_idx  = pz_idx;
            chan_name = 'Pz';
        end
        PSD_WC = zeros(nfft, n_beats);
        f_WC   = zeros(nfft, 1);
        PEEG_raw = zeros(numel(name_bands), n_beats);

        data_norm = zscore(EEG(chan_idx, :));  % z-score Cz and Pz channel


        for k = 1:n_beats
            % Extract EEG segment within cardiac cycle k
            seg_start = posRPeakECG(sel(k));
            seg_end   = posRPeakECG(sel(k) + 1);
            data_tmp  = data_norm(seg_start:seg_end);

            % Weighted covariance PSD
            [~, PSD_WC(:,k), f_WC(:,1)] = ...
                funct_CalcPSDmodified(data_tmp, Bw, window, nfft, fs, 0);

            % Band power per frequency band
            for ibands = 1:numel(name_bands)
                ind = find(f_WC >= ybands(ibands) & f_WC < ybands(ibands+1));
                PEEG_raw(ibands, k) = mean(PSD_WC(ind, k)) * base_banda(ibands);
            end
        end

        % Limit PSD to 30 Hz
        ind_30Hz = find(f_WC <= 30);
        
        % Outlier correction on EEG band power series
        PEEG_corrected = zeros(numel(name_bands), n_beats);
        for ibands = 1:numel(name_bands)
            corrected = removeOutlier(PEEG_raw(ibands, :));

            % pad if edge sample was deleted
            if length(corrected) < n_beats
                corrected(end+1:n_beats) = corrected(end); % if one sample was deleted from the edge, restore it by repeating the last value
            end


            % Check whether any values in corrected are NaN (clustered outliers)
            if any(isnan(corrected))
                nan_idx = isnan(corrected); % Creates a logical vector of the same length as corrected, where true marks NaN positions and false marks valid positions
                corrected(nan_idx) = interp1(find(~nan_idx), corrected(~nan_idx), ...
                    find(nan_idx), 'nearest', 'extrap'); % NaN replacement through nearest neighbour interpolation
            end

            PEEG_corrected(ibands, :) = corrected;
        end
        
        % Outlier correction on full EEG PSD
        PSD_WC_corrected = zeros(length(ind_30Hz), n_beats);
        for ifreq = 1:length(ind_30Hz)
            corrected = removeOutlier(PSD_WC(ind_30Hz(ifreq), :));
            if length(corrected) < n_beats
                corrected(end+1:n_beats) = corrected(end);
            end
            PSD_WC_corrected(ifreq, :) = corrected;
        end

        % Store results under channel name
        results(isubj).(chan_name).PEEG   = PEEG_corrected;
        results(isubj).(chan_name).PSD_WC = PSD_WC_corrected; % [n_freqs_below30 x n_beats] full PSD matrix
        results(isubj).(chan_name).f_WC   = f_WC(ind_30Hz);   % [n_freqs_below30 x 1] frequency vector
    end
end 

%% Save results
save('C:\Users\Isa\OneDrive\Documenten\Unif 2e master\Masterproef_II\data\derivatives\timeSeries\timeseries_results.mat', 'results', 'name_bands', 'ybands');
fprintf('Done. Results saved.\n');
