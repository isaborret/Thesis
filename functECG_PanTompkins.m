%% detecting R peaks using Pan-Tompkins method
% NOTE: changed Fs to fs_ecg
% NOTE: butter filter needs to be skipped, bc signal is already bandpass filtered. 
% Code from Vergara et al. for their brain-heart pipeline, specifically from study: A. Zaccaro, M. G. Perrucci, E. Parrotta, M. Costantini, and F. Ferri, "Brain-heart interactions are modulated across the respiratory cycle via interoceptive attention," NeuroImage, vol. 262, 2022, Art. no. 119548

function [t_PT, signal_PT, posMassimoPT, tMassimoPT, massimoPT] = functECG_PanTompkins(t, signal_raw, filterOrder, Fc_LP_PT, Fc_HP_PT, MinPeakProminence, MinPeakDistanceTh)

    error(nargchk(2,7,nargin));% requires between 2 and 7 inputs; t and signal_raw are mandatory, the rest have defaults
  
    if nargin < 7, MinPeakDistanceTh = 0.4; end %  minimum distance between peaks in seconds (0.4s = 150 bpm maximum heart rate). Prevents double-detection within a single QRS      
    %if nargin < 6, MinPeakProminence = 0.5e-3; end   % gets overruled later     
    if nargin < 5, Fc_HP_PT = 5; end      %5-15 Hz bandpass filter that allows us to identify the QRS complex
    if nargin < 4, Fc_LP_PT = 15; end
    if nargin < 3, filterOrder = 4; end

    %fs = 1/(t(2)-t(1)); % computes sampling rate from time vector
    fs = 256


    %nomalised
    % designfilt requires cutoff frequencies expressed as a fraction of the Nyquist frequency (Fs/2), not in Hz. 
    % So dividing by fs/2 converts Hz to normalised units where 1.0 = Nyquist.
    Fc_LP_PT_norm = Fc_LP_PT / (fs/2);
    Fc_HP_PT_norm = Fc_HP_PT / (fs/2);

    %% Phase 1: Apply a bandpass filter to reduce physiological noise, 
    % baseline, and 50 Hz interference
    
    design_filter_BP_PT = designfilt('bandpassiir','FilterOrder',filterOrder, ...
    'HalfPowerFrequency1',Fc_HP_PT_norm,'HalfPowerFrequency2',Fc_LP_PT_norm,'DesignMethod','butter');
    signal = filtfilt(design_filter_BP_PT,signal_raw);
    
    %% Phase 2: 5-point derivative to obtain the slope of the QRS complex
    % On the five-point derivative: the derivative measures how steeply the signal is changing at each point. 
    % The QRS complex has the steepest slopes in the entire ECG signal — the signal rises and falls very sharply. P and T waves by contrast change slowly. 
    % So taking the derivative effectively converts the QRS from a sharp peak into a very large value, while P and T waves become small values. 
    % The "five-point" part just means it uses 5 surrounding samples to estimate the slope at each point, 
    % which is more accurate and less sensitive to noise than a simple two-point difference.

    % This is then squared in phase 3 to make everything positive and further amplify the QRS relative to everything else, 
    % making it much easier for findpeaks to reliably detect R-peaks.

    %five point derivative
    for i=3:numel(signal)-2
        fivepderiv(i) = (-signal(i+2)+8*signal(i+1)-8*signal(i-1)+signal(i-2))/(12/fs); 
    end
    
    %% Phase 3: square filter is applied
    % squares the derivative to make all values positive and 
    % amplify large slopes (QRS) relative to small ones (P, T waves)
    
    PT_squaring_signal = fivepderiv.^2;     
    
    %% Phase 4: Applying the Moving Average Filter
    % smooths over a window of 256/7 ≈ 37 samples (~144 ms), 
    % integrating energy over the QRS duration to produce broad peaks suitable for detection
    % => are more easily picked up by findpeaks
    signal_PT = smooth(PT_squaring_signal,fs/7);

    % typical QRS complex lasts roughly 80-120 ms, so a 144 ms window is wide enough 
    % to integrate the full energy of one QRS complex into a single smooth peak 
    % without being so wide that it starts merging adjacent beats together
    
    %% Phase 5: selection of appropriate thresholds, even variable ones
    
    t_PT = t(3:end); % is protecting against edge artifacts
    MinPeakDistance = MinPeakDistanceTh*fs; % enforces minimum 400 ms between peaks 

    MinPeakProminence=3*median(signal_PT); % sets adaptive threshold at 3× the median of the smoothed signal. 
                                           % This is robust to amplitude variations across subjects.
    
    [massimoPT,posMassimoPT] = findpeaks(signal_PT, 'MinPeakProminence',MinPeakProminence,'MinPeakDistance',MinPeakDistance); % detects peaks above the prominence threshold with the minimum distance constraint, returning positions and amplitudes
    tMassimoPT = t_PT(posMassimoPT);
end