function NPV = NPV_structured_coupon_3y(settlement_date, exercised_dates, exercised_discounts, prob_digital, coupon_values)
% Computation of the Structured NPV (Coupon + ER + FC)
% 
% INPUT:
% settlement_date:             dates of start contract
% exercise_dates:              col vec of all dates base
% exercised_discounts:         col vec of all discounts related
% prob_digital:                vec prob digital for the trigger (possible) date
% coupon_values:               vec values of the coupons
% 
% OUTPUT:
% NPV:                         Net Present Value at settlement date
% 
% USES:
% function contracts_yf(settlement_date, dates)

    %% Initialization of the parameters
    
    NPV = 0;

    % Find the required interdelta
    [~, ~, yf_30360] = contracts_yf(settlement_date, exercised_dates);
    interdelta = [yf_30360(1); yf_30360(2:end) - yf_30360(1:end-1)];

    %% Computation of the NPV till 3y
    
    NPV = NPV + prob_digital(4) * interdelta(3) * exercised_discounts(3) *coupon_values(2);

    %% Computation of the NPV till 2y, redemption only at 2 year
    
    NPV = NPV + prob_digital(3) * interdelta(2) * exercised_discounts(2) * coupon_values(1);

    %% Computation of the NPV till 2y, redemption only at 1 year
    
    NPV = NPV + prob_digital(2) * interdelta(1) * exercised_discounts(1) * coupon_values(1);

    %% Computation of the NPV till 3y
    
    NPV = NPV + prob_digital(1) * ...
        (interdelta(1) * exercised_discounts(1) *coupon_values(1) + interdelta(2) * exercised_discounts(2) * coupon_values(1));
    
end % function NPV_structured_coupon_3y