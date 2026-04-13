% Takes a data series, optional R-peak positions, and optional interpolation method. 
% Returns the corrected data and optionally corrected R-peak positions
function [correctedData, posRPeaks_corrected] = removeOutlier(data,posRPeaks,method) 
    %fprintf('method = %s, hasPosRPeaks = %d\n', method, hasPosRPeaks);
    if nargin < 3
        method = 'nearest'; % Default interpolation method is nearest-neighbour
    end
    if nargin < 2
        posRPeaks = [];
        hasPosRPeaks = false;
    else
        hasPosRPeaks = true;
    end 
   %%% detect outlier
   mdata = mean(data);
   sdata = std(data);
   posOutlier = find((data > mdata+3*sdata) | (data < mdata-3*sdata)); % Computes mean and standard deviation of the entire series, 
   % then finds indices of all values more than 3 SD from the mean
   if hasPosRPeaks
       posRPeaks_corrected = posRPeaks;
   end
   if ~isempty(posOutlier)
       if posOutlier(1) == 1
            data = data(2:end); % If the very first sample is an outlier, it gets deleted entirely rather than interpolated
            posOutlier = posOutlier(2:end);
            if hasPosRPeaks
                posRPeaks_corrected = posRPeaks_corrected(2:end);
            end
        elseif posOutlier(end) == length(data)
            data = data(1:end-1); % Same logic for the last sample
            posOutlier = posOutlier(1:end-1);
            if hasPosRPeaks
                if ~isempty(posRPeaks)
                    posRPeaks(end-1) = [];
                end
                posRPeaks_corrected = posRPeaks;
            end
       end
   end

   %% INTERPOLATION
   X = 1:length(data);
   cleanX = X; cleanX(posOutlier) = []; % removes the outlier positions from this index vector, leaving only the indices of valid samples
   cleanData = data(cleanX); % extracts the corresponding valid values
   correctedData = interp1(cleanX,cleanData,X,method); % interpolates at all original positions X using only the clean anchor points cleanX and cleanData

   if hasPosRPeaks
        varargout{1} = correctedData; % Returns the corrected series
        varargout{2} = posRPeaks_corrected;
    else
        varargout{1} = correctedData;
    end

end