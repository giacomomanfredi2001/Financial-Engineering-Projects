function NPV_bank = E1_computationNPV_bank(spol, protection, settlementDate, exercised_dates, exercised_discounts)
% Computation of the fixed leg (PARTY A)
% 
%INPUT:
% spol:                       spread over libor
% protection:                 protection from the contract
% settlementDate:             settlement date
% exercised_dates:            dates to exercise the contract
% exercised_discounts:        discounts to exercise the contract

    %% Initialization
    
    ACT360 = 2;

    NPV_bank = 0;

    %% Computation of the interdelta

    yf = [0; yearfrac(settlementDate, exercised_dates, ACT360)];
    interdelta = yf(2:end) - yf(1:end-1);

    %% Computation of the NPV - spol

    NPV_bank = NPV_bank + spol * (interdelta')*exercised_discounts;

    %% Computation of the NPV - LIBOR

    NPV_bank = NPV_bank + 1 - exercised_discounts(end);

    %% Computation of the NPV - Protection

    NPV_bank = NPV_bank + (1-protection)*exercised_discounts(end);
    
end % function E1_computationNPV_bank