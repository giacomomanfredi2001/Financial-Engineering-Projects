function spot_vols = computation_spot_vols(flat_vols, capsPrices, strikes, ...
    exercised_dates, exercised_discounts, settlementDate, number_years, yearlyPayments)
% Computation of the spot volatilities from the Cap prices, the function
% after computing the delta caplet from the prices previously obtained
% tries to calibrate the spot surface volatility for each strike and each
% signed maturity through the function fzero()
% 
%INPUT:
% flat_vols:            flaat volatilities from the market
% capsPrices:           prices of the Cap computed with the flat volatilities
% strikes:              strikes given from the market
% exercised_dates:      dates where to exercise the Euribor 3m
% exercised_discounts:  discounts of the Euribor 3m
% settlementDate:       settlement date of the contract
% number_years:         years where to price the caplet
% yearlyPayments:       number of Euribor in a year
% 
%USES:
% function caps_price(discounts_payments, yf_365, interdelta, Libor_rates, strikes, flat_vol, flag)
    
    %% Conventions
    ACT360 = 2;
    ACT365 = 3;

    %% Quantities of interest
    % Computation of the quantities necessary to get back the spot
    % volatilities from the delta prices of the caplets.

    yearfracs = yearfrac(settlementDate, exercised_dates, ACT360);
    yearfracs_365 = yearfrac(settlementDate, exercised_dates, ACT365);
    interdelta = yearfracs(2:end) - yearfracs(1:end-1);
    forward_discounts = exercised_discounts(2:end) ./ exercised_discounts(1:end-1);
    Libor_rates = (1 ./ forward_discounts - 1) ./ interdelta;

    % Computing delta prices
    delta_prices = capsPrices(2:end,:) - capsPrices(1:end-1,:);
    
    %% Computing spot vols

    % Initialization of the spot volatilities matrix
    spot_vols = zeros(number_years(end)*4 -1,length(strikes));

    % Computation of the 1y
    % The first year owns this structure since it lacks of possibility to
    % be interpolated, the volatility of all the caplets at 1y are then
    % assumed "near" the actual time. Consequently they're assumed equal to
    % the related flat volatilities.
    spot_vols(1:yearlyPayments-1, :) = repmat(flat_vols(1, :), 3, 1);

    % Computation from the 1y to the 10y
    % In order to interpolate the value from 1 to 10y we calibrated the
    % final spot volatility of the period and then interpolated the others.
    % Since it's each year after the other, we can move from a 4 to 4
    % caplet description.

    for j = 1 : length(strikes)
        for i = yearlyPayments : yearlyPayments : yearlyPayments * (number_years(end-2) - 1)
    
            fun = @(z) delta_prices(i/yearlyPayments, j) - caps_price(exercised_discounts(i+1:i+4), yearfracs_365(i:i+3),...
                                                   interdelta(i:i+3), Libor_rates(i:i+3), strikes(j),...
                                                   interp1([exercised_dates(i); exercised_dates(i+4)], [spot_vols(i-1,j); z], exercised_dates(i+1:i+4)), 2);
            spot_vols(i+3,j) = fzero(@(z) fun(z), spot_vols(i-1, j));
            spot_vols(i:i+2,j) = interp1([exercised_dates(i); exercised_dates(i+4)], [spot_vols(i-1,j); spot_vols(i+3,j)], exercised_dates(i+1:i+3)');
        end
    end
    

    % Computation from the 10y to 12y
    % Same as previous but using 8 caplets from the bi annual period

    for j = 1 : length(strikes)
        i = number_years(end-2)*yearlyPayments;
        
        fun = @(z) delta_prices(i/yearlyPayments,j) - caps_price(exercised_discounts(i+1:i+8), yearfracs_365(i:i+7),...
                                               interdelta(i:i+7), Libor_rates(i:i+7), strikes(j),...
                                               interp1([exercised_dates(i); exercised_dates(i+8)], [spot_vols(i-1,j); z], exercised_dates(i+1:i+8)), 2);
        spot_vols(i+7,j) = fzero(@(z) fun(z), spot_vols(i-1, j));
        spot_vols(i:i+6,j) = interp1([exercised_dates(i); exercised_dates(i+8)], [spot_vols(i-1,j); spot_vols(i+7,j)], exercised_dates(i+1:i+7)');
    end

    
    % Computation from the 12y to 15y
    % Same as previous but using 12 caplets from the tri annual period

    for j = 1 : length(strikes)
        i = number_years(end-1)*4;
        fun = @(z) delta_prices(i/4 - 1,j) - caps_price(exercised_discounts(i+1:i+12), yearfracs_365(i:i+11),...
                                               interdelta(i:i+11), Libor_rates(i:i+11), strikes(j),...
                                               interp1([exercised_dates(i); exercised_dates(i+12)], [spot_vols(i-1,j); z], exercised_dates(i+1:i+12)), 2);
        spot_vols(i+11,j) = fzero(fun, spot_vols(i-1, j));
        spot_vols(i:i+10,j) = interp1([exercised_dates(i); exercised_dates(i+12)], [spot_vols(i-1,j); spot_vols(i+11,j)], exercised_dates(i+1:i+11)');
    end
    
end % function computation_spot_vols