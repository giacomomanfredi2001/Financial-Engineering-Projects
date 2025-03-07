function yearlyDiscounts = yearlySwapsDiscounts (dates, discounts, swapsDates)
% Computation of swaps yearly discounts
%
%INPUT
% dates:        dates from bootstrap
% discounts:    discounts from bootstrap
% swapsDates:   dates of swaps payments

    % Conventions
    ACT365 = 3;
    
    % Quantities of interest
    yf = yearfrac(dates(1), swapsDates, ACT365);
    
    % Interpolation of zero rates
    zRates = zeroRates(dates, discounts)/100;
    yearlyRates = interp1(dates(2:end), zRates(2:end), swapsDates);
    
    % Computation of discounts
    yearlyDiscounts = [1; exp(-yf.*yearlyRates)];

end % function yearlySwapsDiscounts