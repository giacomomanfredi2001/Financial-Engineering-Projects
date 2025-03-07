function FV = FV_risky_bond(cf_schedule, Q, ZC_curve, R) 
% Computation of the Forward value with transition matrix
% The method is generalized with as many years coupon as possible but with
% only 1 year jump
% 
%INPUT:
% cf_schedule:        table of cash flows of corp. bonds
% Q:                  transition matrix
% ZC_curve:           table of zero-coupon rates 
% R:                  recovery rate

    %% Computation of the B(1y, tk)

    % Computation of B(0y, tk)
    ZC_rates_0y = interp1(ZC_curve(:, 1), ZC_curve(:, 2), cf_schedule(:, 1), 'linear');
    ZC_Bond_0y = exp(- ZC_rates_0y.*cf_schedule(:, 1));

    % Computation of B(1y, tk)
    ZC_1y = ZC_Bond_0y(find(cf_schedule(:, 1) == 1));
    indexes = find(cf_schedule(:, 1) > 1);

    ZC_Bond_1y = ZC_Bond_0y(indexes)./ZC_1y;

    % Compute survival probability
    surv_probs = zeros(length(indexes) + 1, length(Q) - 1);

    surv_probs(1, :) = 1;
    surv_probs(end, :) = 1 - Q(1:end-1, end);

    % INTENSITY CONSIDER CONSTANT OVER THE ENTIRE TIME PERIOD 
    % NOT PIECEWISE CONSTANT EACH YEAR
    surv_probs(2:end-1, :) = surv_probs(end, :).^(cf_schedule(indexes(1:end-1), 1) - 1);

    % Compute the FV
    FV = zeros(length(Q),1);

    for i = 1:length(Q)-1 
        recovery = 100 * R * ZC_Bond_1y' * (surv_probs(1:end-1, i)-surv_probs(2:end, i));
        coupon = (cf_schedule(indexes, 2).*ZC_Bond_1y)' * surv_probs(2:end, i);

        FV(i) = coupon + recovery;
    end

    % Case default between 0 and 1
    FV(end) = 100 * R;
  
end %function FV_risky_bond