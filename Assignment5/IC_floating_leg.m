function payoff = IC_floating_leg(ES, weight, protection, maturity_discount)
% Computation of the IC of the floating leg
% 
%INPUT:
% ES:                             table of underlying
% weight:                         single weight of the return
% protection:                     protection of the contract
% maturity_discount:              discount at final maturity 

    %% Computation of the underlying value at maturity

    % Introduction
    returns = ES(:, 2:end, :)./ES(:, 1:end-1, :);

    % Underlying
    St = sum(weight * sum(returns, 3), 2);

    %% Computation of the payoff

    payoff = maturity_discount*max(St - protection, 0);

end % function IC_floating_leg