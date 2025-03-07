function [datesSet, ratesSet] = completionSwapsCurve(datesSet, ratesSet, max_year)
% Completation of the swap curve for dates and rates
% 
%INPUT:
% datesSet:         struct with all the dates                 
% ratesSet:         struct with all the rates
% max_year:         max horizon of time
% 
%OUTPUT:
% datesSet:         updates struct with all the swap horizon
% ratesSet:         updates struct with all the swap horizon

    % Calculate the dates yearly
    dates_swaps_complete = zeros(max_year, 1);
    
    for i = 1:length(dates_swaps_complete)
        dates_swaps_complete(i) = datenum(busdate(datetime(datesSet.settlement, 'ConvertFrom', 'datenum') - caldays(1) + calyears(i), 1, eurCalendar));
    end
    
    % Interpolation for the rates
    mid_rates_swaps_complete = interp1(datesSet.swaps, ratesSet.swaps, dates_swaps_complete, 'spline');
    
    % Upgrade of the datesSet and ratesSet
    datesSet.swaps = dates_swaps_complete;
    ratesSet.swaps = mid_rates_swaps_complete;

end % function completionSwapsCurve