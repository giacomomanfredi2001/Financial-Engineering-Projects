function NPV = computeNPVswap(settlementDate, dates, discounts, interp_dates, fixed_rate)
% Computation of the Net Present Value for a swap
% 
%INPUT
% settlementDate:           settlement date of the bootstrap
% dates:                    dates from the bootstrap
% discounts:                discounts from the bootstrap 
% interp_dates:             dates to be interpolated for fixed leg
% fixed_rate:               fixed rate S

    %% Conventions

    ACT30360 = 6;

    %% Computation of the Floating Leg

    B_swap = interpolationDiscounts(dates, discounts, interp_dates);
    
    NPV_float = 1 - B_swap(end);
    
    %% Computation of the Fixed Leg

    delta = [0; yearfrac(settlementDate, interp_dates, ACT30360)];
    interdelta = delta(2:end) - delta(1:end-1);
    
    NPV_fixed = fixed_rate * interdelta'*B_swap;
    
    NPV = NPV_float - NPV_fixed;

end %function computeNPVswap