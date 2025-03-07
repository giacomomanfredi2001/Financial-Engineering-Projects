function prices = prices_matrix(flat_volatilities, strikes, exercise_dates, exercise_discounts, settlement_date, yearlyPayments)
% Computation of the caps' prices from the quotations
% 
%INPUT:
% flat_volatilities:           matrix of quoted flat volatilities     
% strikes:                     column vector of strikes for which we have
%                              the quoations
% exercise_dates:              set of payment dates
% exercise_discounts:          discounts in the payment dates
% settlement_date:             settlement date
% yearlyPayments:              number of caplets per year
%
%OUTPUT:
% prices:                      caps price for the quoted caps
% 
%USES:
% function caps_price(discounts_payments, delta_payments, interdelta, Libor_rates, strikes, flat_vol)

%REMARK: all input vectors must be column vectors

    %% Setting day count convention
    ACT360 = 2;
    ACT365 = 3;
    
    %% Computing quarterly Libor rates

    % Computing the interdelta, forward discounts and libor rates

    yearfracs = yearfrac(settlement_date, exercise_dates, ACT360);
    interdelta = yearfracs(2:end) - yearfracs(1:end-1);
    forward_discounts = exercise_discounts(2:end) ./ exercise_discounts(1:end-1);
    Libor_rates = (1 ./ forward_discounts - 1) ./ interdelta;
    
    yearfracs_365 = yearfrac(settlement_date, exercise_dates, ACT365);

    %% Computation of the caplets

    % Defining a vector that contains the number of caplets that are contained
    % for each cap in our data
    number_of_caplets = [(yearlyPayments - 1 : yearlyPayments : yearlyPayments*10 - 1)'; 47; 59; 79; 99; 119];
    
    % Computing prices
    prices = zeros(size(flat_volatilities));

    for i = 1 : size(flat_volatilities, 1)

        % Computing cap prices vectorially for row(strikes)
        cap = caps_price(exercise_discounts(2:number_of_caplets(i)+1), yearfracs_365(1:number_of_caplets(i)), ...
                         interdelta(1:number_of_caplets(i)),Libor_rates(1:number_of_caplets(i)),...
                         strikes',flat_volatilities(i, :), 1);

        prices(i,:) = cap';
    end

end % function prices_matrix