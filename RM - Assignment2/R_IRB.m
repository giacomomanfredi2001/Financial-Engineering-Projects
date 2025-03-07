function R = R_IRB(PD) 
% Compute the correlations for large corporates
% 
%INPUT:
% PD:         probability of default

    % Computation of the R
    Rmin = 0.12; Rmax = 0.24;
    k = 50;

    R = Rmin*(1-exp(-k*PD))/(1-exp(-k)) + Rmax*exp(-k*PD)/(1-exp(-k));
end % function R_IRB