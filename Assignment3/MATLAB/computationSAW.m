function spreadASW = computationSAW(dates, discounts)
% Computation of the Asset Swap Spread through the theoretical formula of
% the Asset Swap Package
% 
%INPUT
% dates:        vector of dates of the interpolated rates
% discounts:    related vector of discounts of the interp. rates
% 
%USES:
% zeroRates(dates, discounts)

    %% Convention standards
    std_act_360 = 2;
    std_act_365 = 3;
    std_30_360 = 6;

    %% Introduction of the given quantities
    c_bar_0 = 1.01;
    coupon = 0.039;
    
    %% Fixed leg part [Swap 1y; Swap 2y; Swap 3y]
    
    % Delta computation
    yf_fixedleg = [0; 1; yearfrac(dates(1), dates(12:13), std_30_360)];  
    yf_interdelta_fixed = yf_fixedleg(2:end)- yf_fixedleg(1:end-1);
    
    % Computation of the first discount B(t0, ti)
    zero_rates = zeroRates(dates, discounts)/100;
    delta_ACT_365_total = yearfrac(dates(1), dates(2:end), std_act_365);
    zr_1y_swap = interp1(delta_ACT_365_total, zero_rates(2:end), 1, 'linear');
    
    discounts_fixedleg = [exp(-zr_1y_swap*1); discounts(12:13)];
    
    % Computation of C(0)
    face_value = discounts_fixedleg(end);
    c_0 = face_value + coupon * yf_interdelta_fixed' * discounts_fixedleg;
    
    %% Floating leg part [each 3 months]
    
    % Delta computation
    yf_floatingleg = [0; 1; yearfrac(dates(1), dates(12:13), std_act_360)];
    total_len = length(yf_floatingleg);
    
    for i = 2 : total_len
        temp_vec = (linspace(yf_floatingleg(i-1), yf_floatingleg(i), 5))';
        yf_floatingleg = [yf_floatingleg; temp_vec(2:end-1)];
    end
    
    yf_floatingleg = sort(yf_floatingleg);
    yf_interdelta_floating = yf_floatingleg(2:end)- yf_floatingleg(1:end-1);
    
    % Interpolation of the discounts
    zr_floatingleg = interp1(delta_ACT_365_total, zero_rates(2:end), yf_floatingleg, 'linear');
    discounts_floatingleg = [exp(-zr_floatingleg.*yf_floatingleg)];
    
    BPV_floatingleg_0 = yf_interdelta_floating' * discounts_floatingleg(2:end);
    
    
    %% Computation of the Spread over Asset Swap
    
    spreadASW = (c_0 - c_bar_0) / BPV_floatingleg_0;

end % function computationSAW