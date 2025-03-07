function coupon_put = CouponPut(payment_dates, payment_discounts, year, yearly_payments, coupon, ...
    strike_CC, a, sigma)
% Computation of the coupon put through the Jamshidian approach
% 
% INPUT:
% payment_dates:                   exercised dates of the swap
% payment_discounts:               exercised discounts of the swaps
% year:                            year in which the swap option is triggered
% yearly_payments:                 payments of the coupon each year
% coupon:                          value of the coupon
% strike_CC:                       strike of the coupon call
% a:                               speed of mean-reversion of the OU in HW
% sigma:                           volatility of the OU in HW
% 
% OUTPUT:
% coupon_put:                      value of the payer swaption
% 

    %% Initialization of the parameters

    % Reduction of the year fractions required
    yf_365 = payment_dates(year * yearly_payments : end, 2);

    % Computation of the interdelta 
    % (cut those up to 2y since we have a non-call 2y)
    interdelta_30360 = payment_dates(2:end, 4) - payment_dates(1:end-1, 4);
    interdelta_30360 = interdelta_30360(year * yearly_payments : end);

    %% Calibration X_star

    coupon_i = coupon * interdelta_30360;
    coupon_i(end) = coupon_i(end) + 1;

    % Creation of the cutted discounts
    payment_discounts = payment_discounts(year * yearly_payments : end);
    forward_discounts = payment_discounts(2:end)/payment_discounts(1);

    % Function handle for the HJM sigma
    sigma_HJM = @(u, T) sigma/a * (1 - exp(-a*(T-u)));

    % Computation of the sigma_tau
    sigma_tau = sigma_HJM(0, yf_365(2:end) - yf_365(1));

    % Computation of the integral value
    integrand = @(u) sigma_HJM(u, yf_365(2:end)).^2 - sigma_HJM(u, yf_365(1)).^2;
    integrale = integral(@(u) integrand(u), 0, yf_365(1), 'ArrayValued', true);
        
    % Computation of the B(ti, ti, ti + tau)
    triple_discounts = @(x) forward_discounts .* exp(-x .* sigma_tau./sigma - 0.5 .* integrale);

    % Computation of the P
    P = @(x) coupon_i' * triple_discounts(x);

    % Computation of the x_star
    x_star = fzero(@(x) P(x) - strike_CC, 0);

    %% Pricing of the HJM ZCB Call Coupon

    % Computationn of B_star
    B_star = triple_discounts(x_star);

    % Creation of the Black - alike formula
    integrale = integral(@(u) (sigma_HJM(u, yf_365(2:end)) - sigma_HJM(u, yf_365(1))).^2, 0, yf_365(1), 'ArrayValued', true);
    V = sqrt(1/yf_365(1) * (integrale));

    d1 = log(forward_discounts./B_star)./(V .* sqrt(yf_365(1))) + 0.5 .* V * sqrt(yf_365(1));
    d2 = log(forward_discounts./B_star)./(V .* sqrt(yf_365(1))) - 0.5 .* V * sqrt(yf_365(1));

    ZCB_coupon_call = payment_discounts(1) * (forward_discounts.*normcdf(d1) - B_star .* normcdf(d2));

    % Closure of the Jamshidian approach
    coupon_call = coupon_i' * ZCB_coupon_call;

    %% Put/Call parity for the Coupon Bonds

    coupon_put = coupon_call - coupon_i' * payment_discounts(2:end) + payment_discounts(1) * 1;

end % function CouponPut