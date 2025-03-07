function price = European_swaption_price (dates, discounts, bermudanYearFrac, interdeltaSwaps, yearlyDiscounts, a_HW, sigma_HW, N, strike)
% European Swaption Pricing via Hull-White
%
%INPUT
% dates:                    dates vector from bootstrap
% discounts:                discounts vector from bootstrap
% bermudanYearFrac:         vector of Bermudan reset dates [1:10]
% interdeltaSwaps:          vector of fixed leg interdeltas
% yearlyDiscounts:          yearly discounts of swaps
% a_HW:                     Hull-White parameter
% sigma_HW:                 Hull-White sigma
% N:                        number of time steps
% strike:                   Swaption strike
%
%USES
% function sigmaAssociatedHJM (sigmaHW, aHW, u, t)
% function computation_forward_discount (dates, discounts, yfTree)


    %% Tree parameters
    
    % Computing the time grid as a column vector up to year 9 since it is
    % the last date where we can exercise our option
    T = bermudanYearFrac(end-1);
    dt = T/N;
    time_grid = (0 : dt : T)';

    % Computation of model parameters
    mu_hat = 1 - exp(-a_HW*dt);
    sigma_hat = sigma_HW * sqrt( (1-exp(-2*a_HW*dt)) / (2*a_HW) );
    sigma_hat_star = sigma_HW/a_HW * sqrt(dt - 2*(1-exp(-a_HW*dt))/a_HW + (1-exp(-2*a_HW*dt))/(2*a_HW));

    % Computation of delta X Hull-White
    deltaX = sigma_hat*sqrt(3);
    
    % L limits computation
    % From the paper we take l_max as the first integer greater than the
    % lower bound
    l_max = ceil( (1-sqrt(2/3)) / mu_hat );

    % Computation of the vertical grid
    l_grid = (l_max: -1: -l_max)';
    X_grid = l_grid*deltaX;
   
    
    %% Tree probabilities
    
    % Computation of the tree probabilities for each scheme
    % C ---> first line
    % B ---> last line
    % A ---> middle lines
    p_Up = [1/2*(7/3-3*l_grid(1)*mu_hat+(l_grid(1)*mu_hat)^2);...
            1/2*(1/3-l_grid(2:end-1)*mu_hat+(l_grid(2:end-1)*mu_hat).^2);...
            1/2*(1/3+l_grid(end)*mu_hat+(l_grid(end)*mu_hat)^2)];
    
    p_Mid = [-1/3+2*l_grid(1)*mu_hat-(l_grid(1)*mu_hat)^2;...
             2/3-(l_grid(2:end-1)*mu_hat).^2; ...
             -1/3-2*l_grid(end)*mu_hat-(l_grid(end)*mu_hat)^2];
        
    p_Down = [1/2*(1/3-l_grid(1)*mu_hat+(l_grid(1)*mu_hat)^2);...
              1/2*(1/3+l_grid(2:end-1)*mu_hat+(l_grid(2:end-1)*mu_hat).^2);...
              1/2*(7/3+3*l_grid(end)*mu_hat+(l_grid(end)*mu_hat)^2)];
   
          
    %% Computation of the swaption price at the 9th year

    % Computing the forward discount from the bootstrap between the 9th and
    % 10th year 
    initialFwdDiscounts = yearlyDiscounts(end)./yearlyDiscounts(end-1);
    
    % Computation of fwd discounts B(9Y; 9Y, 10Y)
    I = integral(@(u) sigmaAssociatedHJM(sigma_HW,a_HW,u,bermudanYearFrac(end)).^2 - sigmaAssociatedHJM(sigma_HW,a_HW,u,bermudanYearFrac(end-1)).^2, 0, bermudanYearFrac(end-1));
    Discount9Y = initialFwdDiscounts * exp(-X_grid / sigma_HW * sigmaAssociatedHJM(sigma_HW,a_HW,0,bermudanYearFrac(end)-bermudanYearFrac(end-1)) - 1/2 * I);
    
    % Computing Swap rates
    swap_rates = (1-Discount9Y)./(interdeltaSwaps(end)*Discount9Y);

    % Swaption price at 9th year
    vect_old = interdeltaSwaps(end) * Discount9Y .* max(swap_rates - strike, 0);
    
        %% Computing the horizontal part of the tree
    % Initializing the new vector that we will obtain from the continuation
    % value and the early exercise
    vect_new = ones(size(vect_old));

    for i = length(time_grid) - 1 : -1 : l_max + 1

        % Discounts for the time step
      
        % Computing the discount B( t(0); t(i), t(i+1) ) 
        fwdDiscounts_t0 = computation_forward_discount (dates, discounts, time_grid(i:i+1))';

        % Computing the integral in the exponent
        single_step_integral = integral(@(u) sigmaAssociatedHJM(sigma_HW,a_HW,u,time_grid(i+1)).^2 - sigmaAssociatedHJM(sigma_HW,a_HW,u,time_grid(i)).^2, 0, time_grid(i));

        % Computing the discounts
        single_step_discounts = fwdDiscounts_t0* exp(-X_grid*sigmaAssociatedHJM (sigma_HW, a_HW, 0, dt)/sigma_HW - 0.5*single_step_integral);
    
        % Computing stochastic discounts
        stochastic_discount_Up = single_step_discounts .*...
                                 [exp(-0.5 * sigma_hat_star^2 - sigma_hat_star/sigma_hat*( mu_hat*X_grid(1)) );...
                                 exp(-0.5 * sigma_hat_star^2 - sigma_hat_star/sigma_hat*(deltaX + mu_hat*X_grid(2:end-1)));...
                                 exp(-0.5 * sigma_hat_star^2 - sigma_hat_star/sigma_hat*(2*deltaX + mu_hat*X_grid(end)))];
        stochastic_discount_Mid = single_step_discounts .*...
                                  [exp(-0.5 * sigma_hat_star^2 - sigma_hat_star/sigma_hat*( -deltaX + mu_hat*X_grid(1)) );...
                                  exp(-0.5 * sigma_hat_star^2 - sigma_hat_star/sigma_hat*(mu_hat*X_grid(2:end-1)));...
                                  exp(-0.5 * sigma_hat_star^2 - sigma_hat_star/sigma_hat*(deltaX + mu_hat*X_grid(end)))];
        stochastic_discount_Down = single_step_discounts .*...
                                   [exp(-0.5 * sigma_hat_star^2 - sigma_hat_star/sigma_hat*( -2*deltaX + mu_hat*X_grid(1)) );...
                                   exp(-0.5 * sigma_hat_star^2 - sigma_hat_star/sigma_hat*(-deltaX + mu_hat*X_grid(2:end-1)));...
                                   exp(-0.5 * sigma_hat_star^2 - sigma_hat_star/sigma_hat*( mu_hat*X_grid(end)))];
    
        
        % Computing the new continuation value
        vect_new(1) = stochastic_discount_Up(1) * p_Up(1)*vect_old(1) + stochastic_discount_Mid(1)*p_Mid(1)*vect_old(2) + stochastic_discount_Down(1)*p_Down(1)*vect_old(3);
        
        vect_new(2:end-1) = stochastic_discount_Up(2 : end-1) .* p_Up(2:end-1).*vect_old(1:end-2) + stochastic_discount_Mid(2:end-1) .* p_Mid(2:end-1).*vect_old(2:end-1) + stochastic_discount_Down(2:end-1) .* p_Down(2:end-1).*vect_old(3:end);
        
        vect_new(end) = stochastic_discount_Up(end) * p_Up(end) * vect_old(end-2) + stochastic_discount_Mid(end) * p_Mid(end)*vect_old(end-1) + stochastic_discount_Down(end) * p_Down(end) * vect_old(end);
        
        vect_old = vect_new;

    end

    
    %% Computing the triangular part of the tree

    for i = l_max : -1 : 1

        % Computing the discount B( t(i); t(i), t(i+1) ) 
        fwdDiscounts_t0 = computation_forward_discount (dates, discounts, time_grid(i:i+1))';

        % Computing the integral in the exponent
        single_step_integral = integral(@(u) sigmaAssociatedHJM(sigma_HW,a_HW,u,time_grid(i+1)).^2 - sigmaAssociatedHJM(sigma_HW,a_HW,u,time_grid(i)).^2, 0, time_grid(i));

        % Computing the discounts
        single_step_discounts = fwdDiscounts_t0 * exp(-X_grid(-i+l_max+2:end-(l_max-i+1)) * sigmaAssociatedHJM (sigma_HW, a_HW, 0, dt)/sigma_HW - 0.5*single_step_integral);

        % Computing stochastic discounts
        stochastic_discount_Up = single_step_discounts .*  exp(-0.5 * sigma_hat_star^2 - sigma_hat_star/sigma_hat * (deltaX + mu_hat * X_grid(-i+l_max+2 : end-(l_max-i+1))));
        stochastic_discount_Mid = single_step_discounts .* exp(-0.5 * sigma_hat_star^2 - sigma_hat_star/sigma_hat * (mu_hat * X_grid(-i+l_max+2 : end-(l_max-i+1))));
        stochastic_discount_Down = single_step_discounts .* exp(-0.5 * sigma_hat_star^2 - sigma_hat_star/sigma_hat * (-deltaX + mu_hat*X_grid(-i+l_max+2 : end-(l_max-i+1))));

        % Computing the continuation values
        vect_new = stochastic_discount_Up .* p_Up(-i+l_max+2 : end-(l_max-i+1)) .* vect_old(1 : end-2) + stochastic_discount_Mid .* p_Mid(-i+l_max+2 : end-(l_max-i+1)) .* vect_old(2 : end-1) + stochastic_discount_Down .* p_Down(-i+l_max+2 : end-(l_max-i+1)) .* vect_old(3:end);

        vect_old = vect_new;

    end
    price = vect_new;


end % function European_swaption_price 
