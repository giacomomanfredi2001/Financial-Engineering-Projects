function vega= VegaKI(F0,K,KI,B,T,sigma,N,flagNum)
%Computation of the Vega Greek for the hedging
%
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% KI:    barrier
% T:     time-to-maturity
% sigma: volatility
% N:     number of steps
% flagNum: 1 CRR, 2 MC, 3 Exact
%
%USES
% function EuropeanOptionKICRR(F0,K, KI,B,T,sigma,N)

switch (flagNum)
    case 1  % CRR
         
        deltaSigma = 0.05;       % Increments of variance

        % Computation of Centered Euler

        Nsim = 10;
        optionPriceVectorForward = [];
        optionPriceVectorBackward = [];

        for i = 1:Nsim
            optionPriceVectorForward = [optionPriceVectorForward EuropeanOptionKICRR(F0,K, KI,B,T,sigma + deltaSigma,N)];
            optionPriceVectorBackward = [optionPriceVectorBackward EuropeanOptionKICRR(F0,K, KI,B,T,sigma - deltaSigma,N)] ;
    
            optionPriceVectorComplete = (optionPriceVectorForward - optionPriceVectorBackward)/(2*deltaSigma);
        end

        vega = B * mean(optionPriceVectorComplete);

    case 2  % MC
        
        deltaSigma = 0.005;       % Increments of variance
        g = randn(N, 1);          % Normal std distribution

        % Computation of Centered Euler

        optionPriceVectorForward = F0 .* exp(-0.5 * (sigma+deltaSigma)^2 * T - (sigma+deltaSigma) .* sqrt(T) .* g);
        Barrier = optionPriceVectorForward >= KI;
        optionPriceVectorForward = max(optionPriceVectorForward.*Barrier - K, 0);

        optionPriceVectorBackward = F0 .* exp(-0.5 * (sigma-deltaSigma)^2 * T - (sigma-deltaSigma) .* sqrt(T) .* g);
        Barrier = optionPriceVectorBackward >= KI;
        optionPriceVectorBackward = max(optionPriceVectorBackward.*Barrier - K, 0);

        optionPriceVectorComplete = (optionPriceVectorForward - optionPriceVectorBackward)/(2*deltaSigma);

        vega = B * mean(optionPriceVectorComplete);

    case 3  % Closed     
        d2Barrier = (log(F0./KI) - T*sigma^2/2)/(sigma*sqrt(T));

        vegaCall = B.*F0.*exp(-(d2Barrier.^2)/2)/sqrt(2*pi)*sqrt(T);
        vegaDigital = B*(KI-K).*exp(-(d2Barrier.^2)/2)/sqrt(2*pi).*(-sigma^2*T-log(F0./KI)+sigma^2/2*T)/(sigma^2*sqrt(T));
        vega = vegaCall + vegaDigital;
 
    otherwise
end

end %function VegaKI
