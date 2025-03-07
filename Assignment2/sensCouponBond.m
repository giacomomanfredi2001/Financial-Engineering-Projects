function  MacD = sensCouponBond(setDate, couponPaymentDates, fixedRate, dates, discounts)
% Computation of the Macauly Duration for the IB Fixed Coupon
%
%INPUT
% setDate:                   settlement date
% couponPaymentsDates:       payments dates of the coupons
% fixedRate:                 fixed rate of the coupon
% dates:                     dates given by the bootstrap
% discounts:                 discounts given by the bootstrap

    % Computation of the discounts
    index = find(dates == couponPaymentDates(end));
    B_IB_vector = discounts(index-length(couponPaymentDates)+1: index);
    
    % Computation of the interdates deltas
    delta_vector = zeros(length(couponPaymentDates)+1, 1);
    delta_vector(1) = setDate;
    delta_vector(2:end) = couponPaymentDates;
    
    delta_yf_IB = yearfrac(delta_vector(1:end-1), delta_vector(2:end), 6);

    % Computation of the plain delta(t0, ti)
    yf_IB_plain = yearfrac(setDate, couponPaymentDates, 6);

    % Creation of the coupon values
    ci = delta_yf_IB*fixedRate;
    ci(end) = ci(end) + 1;

    % Computation of the Macaulay
    MacD = ((ci.*yf_IB_plain)'*B_IB_vector)/(ci'*B_IB_vector);

end %function sensCouponBond