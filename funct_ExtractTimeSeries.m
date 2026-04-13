%% Extracting the time series
% gives three outputs:
    % series — RR intervals in milliseconds, length N-1
    % tSeries — timestamp of each interval's starting R-peak
    % beats — beat index vector, used mainly for plotting
% Code from Vergara et al. for their brain-heart pipeline, specifically from study: A. Zaccaro, M. G. Perrucci, E. Parrotta, M. Costantini, and F. Ferri, "Brain-heart interactions are modulated across the respiratory cycle via interoceptive attention," NeuroImage, vol. 262, 2022, Art. no. 119548

function [beats, tSeries,  series] = funct_ExtractTimeSeries(tMaximaSignal, npoints) 
   
    series_full = diff(tMaximaSignal); % computes the time differences between consecutive R-peak timestamps
                                       % = RR intervals in seconds
   
    %series=series_full;   %leave in seconds
    series = 1000*series_full; %convert in ms
    tSeries = tMaximaSignal(1:end-1); % timestamp assigned to each RR interval is the time of the first R-peak of that interval (not the midpoint). So if you have N R-peaks you get N-1 RR intervals, each timestamped at its starting R-peak.
    beats=[1:1:numel(series_full)]; % counts the beats
   
    if (npoints>0) % Extraction_WC does this is that it wants all detected RR intervals first, and then selects 400 valid ones using the ±25% mean threshold in the while loop. If you called it with npoints = 300 instead, it would just blindly take the first 300 RR intervals regardless of whether they're valid beats, which is not what you want.
        if (numel(series) > npoints)
            series = series(1:npoints);
            beats = beats(1:npoints);
            tSeries=tSeries(1:npoints);
        else
            disp('Less than 300 points!');    
        end   
    end
end