function dirtyValue = PV_risky_bond_Z(z, cf_schedule, ZC_curve)
% Computation of the dirty value of the Z-spread
% 
%INPUT:
% z:                 z-spread at the given year
% cf_schedule:       matrix of [dates coupon_values]
% ZC_curve:          table of zero-coupon rates

    %% Interpolation of the discount values B(0, ti)

    ZC_rates = spline(ZC_curve(:, 1), ZC_curve(:, 2), cf_schedule(:, 1));
    ZC_Bond = exp(- (ZC_rates + z).*cf_schedule(:, 1));

    %% Computation of the Dirty Value - No default part

    dirtyValue = cf_schedule(:, 2)'*ZC_Bond;

end % function PV_risky_bond_Z