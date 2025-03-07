function laplace_trasform = LaplaceTransform (deltaT, k, alpha, w, sigma, flag)
% Computation of the Laplace Transform for a general w: L(w)
%
%INPUT:
% deltaT:                    interval of time
% k:                         volvol of the dynamic
% alpha:                     parameter to define a particular model
% w:                         integrand function
% sigma:                     volatility of the dynamic
% flag:                      [1: NIG, 2:VG]

    if flag == 1
        laplace_trasform = exp(deltaT/k * (1-alpha)/alpha * (1 - (1 + (w.*k.*sigma.^2)/(1-alpha)).^alpha));
    else
        laplace_trasform = exp(-deltaT/k * log(1 + k.*w.*sigma.^2)); 
    end

end % function LaplaceTransform