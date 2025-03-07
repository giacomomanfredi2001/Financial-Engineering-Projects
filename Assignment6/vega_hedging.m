function [CAP_initial_price,x_CAP] = vega_hedging(flat_vols,initial_spot_vols,strikes, datesSet, exercised_dates, ...
    exercised_discounts, yearlyPayments, strike_CAP, CAP_TTM, upfront)
% Computation of the vega hedging through the caps prices
% 
%INPUT:
% flat_vols:             matrix of flat volatilities
% initial_spot_vols:     matrix of initial stop vols
% strikes:               vector of strikes from mkt
% datesSet:              initial struct of dates
% exercised_dates:       dates of the euribor
% exercised_discounts:   discounts of the euribor
% yearlyPayments:        payments per year of the euribor
% strike_CAP:            strike of the cap
% CAP_TTM:               time to maturity of the cap
% upfront:               upfront of the certificate
% 
%OUTPUT:
% CAP_initial_price:     initial price of the CAP before hedging
% x_CAP:                 qty to by of the cap
% 
%USES:
% function caps_price(discounts_payments, yf_365, interdelta, Libor_rates, strikes, flat_vol, flag)
% function prices_matrix(flat_volatilities, strikes, exercise_dates, exercise_discounts, settlement_date, yearlyPayments)
% function computation_spot_vols(flat_vols, capsPrices, strikes, exercised_dates, exercised_discounts, settlementDate, number_years, yearlyPayments)
% function computationNPV_bank(spol, settlementDate, exercised_dates, exercised_discounts)
% function computationNPVB_wo_upfront(datesSet, exercised_dates, exercised_discounts, spot_vols, strikes, yearlyPayments)

    %% Quantities of interest 
    ACT360 = 2;
    ACT365 = 3;

    yearfracs = yearfrac(datesSet.settlement, exercised_dates, ACT360);
    yearfracs_365 = yearfrac(datesSet.settlement, exercised_dates, ACT365);
    interdelta = yearfracs(2:end) - yearfracs(1:end-1);
    forward_discounts = exercised_discounts(2:end) ./ exercised_discounts(1:end-1);
    Libor_rates = (1 ./ forward_discounts - 1) ./ interdelta;
    
    %% Computing the initial value of the derivative portfolio - cap initial price
    % Computing the spot vols for the strike of the ATM 5Y CAP
    interpolated_spot_vols = spline(strikes, initial_spot_vols(1 : yearlyPayments * CAP_TTM - 1, :),...
                                    strike_CAP);
    
    % Computing the initial price of the ATM 5Y CAP
    CAP_initial_price = caps_price(exercised_discounts(2 : yearlyPayments * CAP_TTM),...
                                   yearfracs_365(1 : yearlyPayments * CAP_TTM - 1),...
                                   interdelta(1 : yearlyPayments * CAP_TTM - 1),...
                                   Libor_rates(1 : yearlyPayments * CAP_TTM - 1),...
                                   strike_CAP, interpolated_spot_vols, 2);
    
    % First of all we need to hedge our certificate from movements of the flat
    % volatilities with a ATM 5Y CAP
    
    %% Bumping volatilities
    % Bumping flat volatilities (total Vega)
    flat_vols_bumped = flat_vols + 1e-4;

    % Computing the new caps' price after the bump
    prices_bumped = prices_matrix(flat_vols_bumped, strikes, exercised_dates, ...
                                  exercised_discounts, datesSet.settlement, yearlyPayments);

    % Computing bumped spot vols
    number_years = [1:10, 12, 15]';
    spot_vols_bumped = computation_spot_vols(flat_vols_bumped, prices_bumped, strikes, ...
                                             exercised_dates, exercised_discounts, datesSet.settlement,...
                                             number_years, yearlyPayments);

    % Computing the interpolated spot vol after the bump
    interpolated_spot_vols_bumped = spline(strikes, spot_vols_bumped(1 : yearlyPayments * CAP_TTM - 1, :),...
                                    strike_CAP);
    
    %% Computing the new NPVs of the portfolio after the bump

    spol = 2 * 1e-2;

    NPV_bank_bumped = computationNPV_bank(spol, datesSet.settlement, exercised_dates(1:15*yearlyPayments),...
                                          exercised_discounts(1:15*yearlyPayments));

    NPV_IB_bumped = upfront + computationNPVB_wo_upfront(datesSet, exercised_dates,...
                                                         exercised_discounts, spot_vols_bumped,...
                                                         strikes, yearlyPayments);

    CAP_bumped = caps_price(exercised_discounts(2 : yearlyPayments * CAP_TTM),...
                            yearfracs_365(1 : yearlyPayments * CAP_TTM - 1),...
                            interdelta(1 : yearlyPayments * CAP_TTM - 1),...
                            Libor_rates(1 : yearlyPayments * CAP_TTM - 1),...
                            strike_CAP, interpolated_spot_vols_bumped, 2);
    
    % Computing the notional of the ATM 5Y CAP imposing the new NPVC equal to
    % the initial NPV of the portfolio
    x_CAP = (NPV_bank_bumped-NPV_IB_bumped) / (CAP_initial_price-CAP_bumped);

end % function vega_hedging