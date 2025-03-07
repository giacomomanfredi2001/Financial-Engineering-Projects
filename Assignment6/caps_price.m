function caps_price = caps_price(discounts_payments, yf_365, interdelta, Libor_rates, strikes, flat_vol, flag)
% Computation of the price for caps that have a given maturity
%
%INPUT:
% discounts_payments:          discounts at payment dates
% yf_365:                      maturity of the payment
% interdelta:                  yearfrac between two following payment dates
% Libor_rates:                 Libor rates between payment dates
% strikes:                     column vector of strikes for which we have
%                              quotes
% flat_vol:                    quoted flat volatilities (row vector if
%                              computed for all the strikes flat vol; column
%                              vector for a single strike spot vol)
% flag:                        1 sigma flat, 2 sigma spot
%
%OUTPUT:
% caps_price:                  column vector of prices 
% 

switch flag
    case 1
        % In this case we compute for all the strikes, given only a single
        % flat vol for each strike, the yearth cap.
        %

        %% Computing the necessary parameters
    
        % Computing the prefactor in Bachelier's formula
        prefactor = discounts_payments .* interdelta;
        
        % Computing d_n
        num_d_n = Libor_rates - strikes;
        den_d_n = repmat(flat_vol, length(Libor_rates), 1) .* sqrt(yf_365);
        d_n = num_d_n ./ den_d_n;
        
        %% Computing the caplet matrix
    
        % In caplet_matrix each the rows are related to the strikes
        caps_price = (prefactor' * (num_d_n.*normcdf(d_n) + den_d_n.*normpdf(d_n)))';

    case 2
        % In this case we have only a single strike but multiple
        % volatilities for the same strike

        % Computing the prefactor in Bachelier's formula
        prefactor = discounts_payments .* interdelta;
        
        % Computing d_n
        num_d_n = Libor_rates - strikes;
        den_d_n = flat_vol .* sqrt(yf_365);
        d_n = num_d_n ./ den_d_n;
        
        %% Computing the caplet matrix
    
        % In caplet_matrix each the rows are related to the strikes
        caps_price = (prefactor' * (num_d_n.*normcdf(d_n) + den_d_n.*normpdf(d_n)));

    case 3
        % In this case we have a vector of strikes and multiple
        % volatilities

        % Computing the prefactor in Bachelier's formula
        prefactor = discounts_payments .* interdelta;
        
        % Computing d_n
        num_d_n = Libor_rates - strikes;
        den_d_n = flat_vol .* sqrt(yf_365);
        d_n = num_d_n ./ den_d_n;
        
        %% Computing the caplet matrix
    
        % In caplet_matrix each the rows are related to the strikes
        caps_price = (prefactor .* (num_d_n.*normcdf(d_n) + den_d_n.*normpdf(d_n)));
end
    
end % function caps_price