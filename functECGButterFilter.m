%% Butterworth filter fro ECG
% removes low-frequency baseline wander (below 0.1 Hz, e.g. breathing artifacts, slow drift) 
% and high-frequency noise (above 20 Hz, e.g. muscle noise) from the ECG, 
% keeping only the frequency range where the QRS complex lives
% Code from Vergara et al. for their brain-heart pipeline, specifically from study: A. Zaccaro, M. G. Perrucci, E. Parrotta, M. Costantini, and F. Ferri, “Brain-heart interactions are modulated across the respiratory cycle via interoceptive attention,” NeuroImage, vol. 262, 2022, Art. no. 119548


function [t, signalFilt] = functECGButterFilter(t, signal, filterOrder,Fc_LP,Fc_HP)

    error(nargchk(2,5,nargin)); %checks that between 2 and 5 input arguments were provided, throws an error otherwise.
    % following lines  set default values if optional arguments aren't provided
    if nargin < 5, Fc_HP = 0.1; end     % highpass at 0.1 Hz
    if nargin < 4, Fc_LP = 20; end      % lowpass at 20 Hz
    if nargin < 3, filterOrder = 4; end % 4th order filter

    Fs = 1/(t(2)-t(1)); % computes sampling frequency from the time vector by taking the inverse of the time step between samples.

    %following normalizes the cutoff frequencies to the Nyquist frequency
    Fc_LP_ECG_norm = Fc_LP / (Fs/2);
    Fc_HP_ECG_norm = Fc_HP / (Fs/2);

    design_filter_BP_ECG = designfilt('bandpassiir','FilterOrder',filterOrder, ...
    'HalfPowerFrequency1',Fc_HP_ECG_norm,'HalfPowerFrequency2',Fc_LP_ECG_norm,'DesignMethod','butter'); % designs a Butterworth bandpass IIR filter with the specified order and cutoff frequencies

    signalFilt = filtfilt(design_filter_BP_ECG,signal); %applies the filter twice, once forward and once backward
end