function rates = interpolationRates(dates, discounts, interp_dates)
% Computation of the interpolated interest rates
% 
% INPUT:
% dates:               dates bootstrapped
% discounts:           discounts bootstrappes
% interp_dates:        dates to interpolate
% 
% OUTPUT:
% rates:               rates interpolated
% 
% USES:
% function zeroRates()

    zRates = zeroRates(dates, discounts)/100;
    rates = interp1(dates, zRates, interp_dates);

end % function interpolationRates