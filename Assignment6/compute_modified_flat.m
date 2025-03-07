function [flat_vol_5y, flat_vol_15y] = compute_modified_flat(flat_vols, datesSet, bps_shift)
% Computing the flat volatilities matrix course-grained
%
%INPUT:
% flat_vols:                 flat volatilities matrix
% datesSet:                  initial struct of the dates
% bps_shift:                 shift of the rates to perform
% 
%OUTPUT:
% flat_vol_5y:               matrix of the flat volatility at 5y
% flat_vol_15y:              value of the flat volatility at 15 y

    %% Quantities of interest
    quoted_years = [1:10 12 15]';

    flat_vol_5y = flat_vols;
    flat_vol_15y = flat_vols;

    %% COarse Grained 5y-15y
    
    % Computing weights for 15y
    weights = interp1([datesSet.swaps(5), datesSet.swaps(15)], [0, 1], datesSet.swaps(quoted_years(5:end)));
    weights_matrix = repmat(weights, 1, 13);

    % Weighted shifting of volatilities
    flat_vol_15y(5:12, :) = flat_vol_15y(5:12, :) + weights_matrix * bps_shift;
    flat_vol_15y(13:end, :) = flat_vol_15y(13:end, :) + bps_shift;
    
    %% Coarse Grained 0-5y
    
    % Computing weights for 5y
    weights = interp1([datesSet.swaps(5), datesSet.swaps(15)], [1, 0], datesSet.swaps(quoted_years(5:end)));
    weights_matrix = repmat(weights, 1, 13);

    % Weighted shifting of volatilities
    flat_vol_5y(5:12, :) = flat_vol_5y(5:12, :) + weights_matrix * bps_shift;
    flat_vol_5y(1:4, :) = flat_vol_5y(1:4, :) + bps_shift;
    
end % function compute_modified_flat