function zRates = zeroRates(dates, discounts)
% Computation of the Zero rate curve
% 
%INPUT
% dates:        vector of all possible dates
% discounts:    vector of all the discounts related to the date

    % Computation of the year fraction
    yf_dates = yearfrac(dates(1), dates, 3);

    % Computation of the zero rates in percentage units
    zRates = - log(discounts)./yf_dates *100;

end %function zeroRates