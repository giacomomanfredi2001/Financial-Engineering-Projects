function NPV_contingent = computation_FtD_contingentleg(dates, discounts, tau_total, indexes_UCG_default, indexes_ISP_default, recovery)
% Computation of the contingent leg for every simulation
% 
%INPUT:
% dates:                vector of dates of the ZR curve
% discounts:            vector of discounts of the ZR curve
% tau_total:            matrix [ISP UCG] containing the defaults time
% indexes_UCG_default:  indexes of the UCG that defaults
% indexes_ISP_default:  indexes of the ISP that defaults
% recovery:             vector of the recoveries of ISP/UCG
% 
%USES:
% function zeroRates(dates, discounts)

    %% Standard creation
    std_act_365 = 3;

    % Zero rates due for the interpolation
    delta_ACT_365_total = yearfrac(dates(1), dates(2:end), std_act_365);
    zero_rates = zeroRates(dates, discounts)/100;

    %% Computation NPV - ISP default
    NPV_contingent_ISP = 0;

    for i = 1: length(indexes_ISP_default)
        zr_interpolated = interp1(delta_ACT_365_total, zero_rates(2:end), tau_total(indexes_ISP_default(i), 1), 'linear');
        discount_default = exp(-zr_interpolated * tau_total(indexes_ISP_default(i), 1));

        if (isnan(discount_default))
            discount_default = 1;
        end

        NPV_contingent_ISP = NPV_contingent_ISP + (1 - recovery(1))*discount_default;
    end

    %% Computation NPV - UCG default
    NPV_contingent_UCG = 0;

    for i = 1: length(indexes_UCG_default)
        zr_interpolated = interp1(delta_ACT_365_total, zero_rates(2:end), tau_total(indexes_UCG_default(i), 2), 'linear');
        discount_default = exp(-zr_interpolated * tau_total(indexes_UCG_default(i), 2));

        if (isnan(discount_default))
            discount_default = 1;
        end

        NPV_contingent_UCG = NPV_contingent_UCG + (1 - recovery(2))*discount_default;
    end

    %% Computation of the total NPV contingent

    NPV_contingent = NPV_contingent_UCG + NPV_contingent_ISP;

end % function computation_FtD_contingentleg