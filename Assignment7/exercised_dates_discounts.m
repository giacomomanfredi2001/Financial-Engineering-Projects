function [payment_dates, reset_dates, payment_discounts, reset_discounts] = exercised_dates_discounts(datesSet, dates, discounts, n_years, yearly_payments)
% Computation of reset dates, payment dates, reset discounts and payment
% discounts
% 
% INPUT:
% n_years:                 number of years for the computation
% yearly_payments:         number of payments for each year (euribor based)
% 
% OUTPUT:
% payment_dates:           dates for the euribor each yeary_payments months
% reset_dates:             annual reset dates (2 days before the related payment)
% payment_discounts:       discounts at the payment dates
% reset_discounts:         discounts at the reset dates
% 
% USES:
% function contracts_yf(settlement_date, dates)
% function interpolationDiscounts(dates, discounts, interp_dates)

    %% Computation of the dates required [RAW NUMBER, YF_365]
    payment_dates = zeros(yearly_payments * n_years, 3);
    reset_dates = zeros(2, 2);
    
    for i = 1:length(payment_dates)    
        payment_dates(i, 1) = datenum(busdate(datetime(datesSet.settlement, 'ConvertFrom', 'datenum') - ...
            caldays(1) + calmonths(12/yearly_payments * i), 1, eurCalendar));
        [payment_dates(i, 3), payment_dates(i, 2), payment_dates(i, 4)] = contracts_yf(datesSet.settlement, payment_dates(i, 1));
    end
    
    for i = 1:length(reset_dates)    
    
        % Reset dates are the dates of payment minor 2 days
        reset_dates(i) = datenum(busdate(datetime(datesSet.settlement, 'ConvertFrom', 'datenum') - ...
            caldays(1 + 2) + calyears(i), 1, eurCalendar));
        [~, reset_dates(i, 2), ~] = contracts_yf(datesSet.settlement, reset_dates(i, 1));
    end
    
    %% Computation of the discounts

    reset_discounts = interpolationDiscounts(dates, discounts, reset_dates(:, 1));
    payment_discounts = interpolationDiscounts(dates, discounts, payment_dates(:, 1));

end % function exercised_dates_discounts