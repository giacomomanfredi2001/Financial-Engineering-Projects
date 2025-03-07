function NPV_fee = computation_FtD_feeleg(spread, interdelta, discounts, tau_total, indexes_UCG_default, indexes_ISP_default, indexes_no_default, yf_swap)
% Computation of the fee leg for every simulation
% 
%INPUT:
% spread:               spread payed at the resets dates
% interdelta:           delta(ti, ti+1) - amount of time between reset dates
% discounts:            B(t0, ti+1) - discounts at the reset dates
% yf_swap:              year fraction of the reset dates
% tau_total:            matrix [ISP UCG] containing the defaults time
% indexes_UCG_default:  indexes of the UCG that defaults
% indexes_ISP_default:  indexes of the ISP that defaults
% indexes_no_default:   indexes of the UCG/ISP that don't default

    %% Computation NPV - No default
    NPV_fee_nodefault = length(indexes_no_default) * spread * interdelta'*discounts;

    %% Computation NPV - Default ISP
    NPV_fee_default_ISP = 0;
    
    for i = 1: length(indexes_ISP_default)
        indicator = zeros(length(discounts), 1);
        tau = tau_total(indexes_ISP_default(i), 1);

        index = find(yf_swap > tau, 1);
        indicator( 1 : index ) = 1;

        NPV_fee_default_ISP = NPV_fee_default_ISP + spread * (indicator .* interdelta)'*discounts;
    end

    %% Computation NPV - Default UCG
    NPV_fee_default_UCG = 0;
    
    for i = 1: length(indexes_UCG_default)
        indicator = zeros(length(discounts), 1);
        tau = tau_total(indexes_UCG_default(i), 2);

        index = find(yf_swap > tau, 1);
        indicator( 1 : index ) = 1;

        NPV_fee_default_UCG = NPV_fee_default_UCG + spread * (indicator .* interdelta)'*discounts;
    end

    %% Final computation of the NPV
    NPV_fee = NPV_fee_default_UCG + NPV_fee_nodefault + NPV_fee_default_ISP;

end % function computation_FtD_feeleg