%% find true R peak
%  function refines each candidate peak by searching for the true maximum
% within a small window around each Pan-Tompkins estimate
% Code from Vergara et al. for their brain-heart pipeline, specifically from study: A. Zaccaro, M. G. Perrucci, E. Parrotta, M. Costantini, and F. Ferri, "Brain-heart interactions are modulated across the respiratory cycle via interoceptive attention," NeuroImage, vol. 262, 2022, Art. no. 119548

function [posRPeakECG, tPosRPeakECG, massimoRPeakECG] = functECG_findRPeak(t, signal, PosMaxECG_PT)

    for i=1:numel(PosMaxECG_PT) % loops over each candidate R-peak position from Pan-Tompkins
        centralPoint = PosMaxECG_PT(i); % takes the candidate position as the centre of the search window
        winwidth = 40; % ±40 sample search window around the candidate
        startPoint = centralPoint-winwidth; % defines window edges, with boundary checks to avoid going outside the signal.
        if (startPoint < 1)
            startPoint =1;
        end
        endPoint = centralPoint+winwidth;
        if (endPoint > numel(signal))
            endPoint = numel(signal);
        end

        arrayTemp = signal(startPoint:endPoint); % extracts the signal segment within the window
        [massimoRtmp,posMassimoRtmp] = max(arrayTemp); % finds the maximum value and its position within the local window
        % = refined R-peak

        massimoRPeakECG(i,1) = massimoRtmp; % stores the amplitude of the refined R-peak
        posRPeakECG(i,1) = posMassimoRtmp+startPoint-1; % the index correction works like this: imagine startPoint = 100 and max finds the peak at local position 15 within arrayTemp. 
        % If you just used 15, that would be meaningless globally. Adding startPoint gives 115, but this overcounts by 1 because position 100 is already local position 1, not local position 0. So you subtract 1: 15 + 100 - 1 = 114. This gives the correct global sample index.

        tPosRPeakECG(i,1) = t(posRPeakECG(i)); % simply looks up the timestamp corresponding to that global sample index in your time vector t. So if sample 114 corresponds to t=0.445 seconds, that becomes the R-peak timestamp.
    end

end