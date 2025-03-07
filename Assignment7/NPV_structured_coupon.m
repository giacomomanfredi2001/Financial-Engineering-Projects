function NPV = NPV_structured_coupon(settlement_date, exercised_dates, exercised_discounts, prob_digital, coupon_values)
% Computation of the Structured NPV (Coupon + ER + FC)
% 
% INPUT:
% settlement_date:             dates of start contract
% exercise_dates:              col vec of all dates base
% exercised_discounts:         col vec of all discounts related
% prob_digital:                vec prob digital for the trigger (possible) date
% coupon_values:               value of the coupons
% 
% OUTPUT:
% NPV:                         Net Present Value at settlement date
% 
% USES:
% function contracts_yf(settlement_date, dates)

    %% Initialization of the parameters
    
    % Find the required interdelta
    [~, ~, yf_30360] = contracts_yf(settlement_date, exercised_dates);
    interdelta = [yf_30360(1); yf_30360(2:end) - yf_30360(1:end-1)];

    %% Computation of the components

    coupon_total = interdelta .* exercised_discounts .* coupon_values;
    
    %% Computation of the NPV

    NPV = prob_digital' * coupon_total;
    
end % function NPV_structured_coupon