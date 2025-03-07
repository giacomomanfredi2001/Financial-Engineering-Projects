function sigmaHJM = sigmaAssociatedHJM (sigmaHW, aHW, u, t)
% Associated Gaussian HJM volatility computation
%
%INPUT
% sigmaHW:        Hull-White model volatility
% aHW:            Hull-White model parameter
% t:              time for the computation
% u:              t0

    sigmaHJM = sigmaHW/aHW * (1-exp(-aHW*(t-u)));

end % function sigmaAssociatedHJM