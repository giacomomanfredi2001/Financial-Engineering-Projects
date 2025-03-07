function NPV = computeNPVswap(settlementDate, datesSet, yf_dates, discounts, fixed_rate)
% Computation of the Net Present Value for a swap
% 
%INPUT
% settlementDate:           settlement date of the bootstrap
% datesSet:                 dates for the fixed leg
% yf_dates:                 year fraction of all the general dates
% discounts:                discounts of all the dates bootstrapped
% fixed_rate:               fixed rate S

    % Computation of the Floating Leg
    yf_swap_final = yearfrac(settlementDate, datesSet(end), 3);
    
    index = find(yf_dates == yf_swap_final);
    B_swap = discounts(index-length(datesSet)+1: index);
    
    NPV_float = 1 - B_swap(end);
    
    % Computation of the Fixed Leg
    delta_vector = zeros(length(datesSet)+1, 1);
    delta_vector(1) = settlementDate;
    delta_vector(2:end) = datesSet;
    
    delta_yf_swap = yearfrac(delta_vector(1:end-1), delta_vector(2:end), 6);
    
    NPV_fixed = fixed_rate* delta_yf_swap'*B_swap;
    
    NPV = NPV_float - NPV_fixed;

end %function computeNPVswap