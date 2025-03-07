function NPV = NPV_swap_digital(settlement_date, exercised_dates, exercised_discounts, ...
    prob_digital, spol, yearly_payments, trigger_year)
% Computation of the Swap NPV (Euribor + spol)
% 
% INPUT:
% settlement_date:             dates of start contract
% exercise_dates:              col vec of all dates base
% exercised_discounts:         col vec of all discounts related
% prob_digital:                vec prob digital for the cases [p1, p1p2,
%                              1-p1p2]
% spol:                        spread over libor of Euribor "i"m
% yearly_payments:             # of payments of the euribor in a year
% trigger_year:                year when the trigger is pulled
% 
% OUTPUT:
% NPV:                         Net Present Value at settlement date
% 
% USES:
% function contracts_yf(settlement_date, dates)

    %% Initialization of the parameters

    % Find the dates to trigger (possible)
    trigger_index = [trigger_year * yearly_payments; length(exercised_discounts)];
    
    % Find the required interdelta
    [yf_360, ~, ~] = contracts_yf(settlement_date, exercised_dates);
    interdelta = [yf_360(1); yf_360(2:end) - yf_360(1:end-1)];
    
    %% Compute the parts to get the NPV
    libor_part = 1 - exercised_discounts(trigger_index);

    spol_part = zeros(length(trigger_index), 1);

    for i = 1:length(trigger_index)
        spol_part(i) = interdelta(1:trigger_index(i))' * exercised_discounts(1:trigger_index(i)) * spol;
    end

    %% Compute the NPV

    NPV = prob_digital' * (libor_part + spol_part);
    
end % function NPV_swap_digital