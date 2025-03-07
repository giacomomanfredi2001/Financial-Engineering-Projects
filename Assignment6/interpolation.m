function B_temp = interpolation(B1, B2, yf1, yf2, yf)
% Interpolation function
% 
%INPUT
% B1:    discount initial point interval
% B2:    discount final point interval
% yf1:   yf initial point interval
% yf2:   yf final point interval
% yf:    yf interpolated point interval

    % Compute the required quantities
    yt1 = - log(B1)/yf1;
    yt2 = - log(B2)/yf2;

    dti = yf2 - yf1;
    dt = yf - yf1;

    % Interpolate the rates
    y = yt1 + (yt2-yt1)*dt/dti;
    B_temp = exp(-yf*y);

end %function interpolation