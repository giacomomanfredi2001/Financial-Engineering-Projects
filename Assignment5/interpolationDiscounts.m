function discounts_interpolated = interpolationDiscounts(dates, discounts, interp_dates)
% Computation of the interpolated discounts
% 
%INPUT:
% dates:               dates of the bootstrap
% discounts:           discounts of the bootstrap
% interp_dates:        dates to interpolate
% 
%USES:
% function zeroRates(dates, discounts)

    %% Conventions

    conv_ACT365 = 3;

    %% Computation of the discounts

    % Computation of the zero rates
    rates = zeroRates(dates, discounts)/100;
    rates_MD = interp1(dates(2:end), rates(2:end), interp_dates);

    yf_MD = yearfrac(dates(1), interp_dates, conv_ACT365);

    discounts_interpolated = [exp(-rates_MD .* yf_MD)];

end % function interpolationDiscounts