function NPV_inv_bank = E1_computationNPV_investment_bank(upfront, ES, weight, alpha, protection, maturity_discount)
% Computation of the floating leg (PARTY B)
% 
%INPUT:
% upfront:                        initial payment
% ES:                             table of underlying
% weight:                         single weight of the return
% alpha:                          partecipation coefficient
% protection:                     protection of the contract
% maturity_discount:              discount at final maturity 

    %% Computation of the underlying value at maturity

    returns = ES(:, 2:end, :)./ES(:, 1:end-1, :);

    % Underlying
    St = sum(weight * sum(returns, 3), 2);

    %% Computation of the payoff

    payoff = alpha * max(St - protection, 0);

    %% Computation of the NPV

    NPV_inv_bank = sum(upfront + maturity_discount*payoff);

end % function E1_computationNPV_investment_bank