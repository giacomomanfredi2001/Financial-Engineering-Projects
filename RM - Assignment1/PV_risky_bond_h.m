function dirtyValue = PV_risky_bond_h(cf_schedule, h_curve, ZC_curve,R)
% Computation of the dirty value with default risk
% 
%INPUT:
% cf_schedule:       matrix of [dates coupon_values]
% h_curve:           possible hazard values
% ZC_curve:          table of zero-coupon rates
% R:                 recovery value

    %% Modification of the hazard rates

    h_curve = h_curve(:, 2);

    %% Interpolation of the discount values B(0, ti)

    ZC_rates = spline(ZC_curve(:, 1), ZC_curve(:, 2), cf_schedule(:, 1));
    ZC_Bond = exp(- ZC_rates.*cf_schedule(:, 1));

    %% Computation of the Dirty Value - No default part

    survival_prob = zeros(1 + length(ZC_Bond), 1);    
    survival_prob(1) = 1;

    times = [0; cf_schedule(:, 1)];
    
    for i = 2 : length(ZC_Bond) + 1
        survival_prob(i) = survival_prob(i-1)*exp(-h_curve( ceil( times(i) ) )* ( times(i) - times(i-1) ) );
    end

    survival_prob = survival_prob(2:end);
    dirtyValue_nodefault = cf_schedule(:, 2)'*(ZC_Bond.*survival_prob);

    %% Computation of the Dirty Value - Default part

    surv_temp = [1; survival_prob];

    dirtyValue_default = 100*R*(surv_temp(1:end-1) - surv_temp(2:end))'*ZC_Bond;

    %% Assembling of the total value

    dirtyValue = dirtyValue_nodefault + dirtyValue_default;

end %function PV_risky_bond_h