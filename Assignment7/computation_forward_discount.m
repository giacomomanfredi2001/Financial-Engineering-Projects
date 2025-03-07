function fwdDiscounts = computation_forward_discount (dates, discounts, yfTree)
% Forward Discounts for each time step of the tree
%
%INPUT
% dates:            dates from bootstrap
% discounts:        discounts from bootstrap
% yfTree:           year fractions of the tree time steps 

    % Conventions   
    ACT365 = 3;
    
    % Year fraction of input dates
    yfDates = yearfrac(dates(1), dates, ACT365);

    % Interpolation of zero rates
    zRates = zeroRates(dates, discounts)/100;
    zRates(1) = zRates(2);
    treeRates = interp1(yfDates, zRates, yfTree);
    
    % Computation of discounts
    treeDiscounts =  exp(-yfTree.*treeRates);
    
    % Computation of fwd discounts
    fwdDiscounts = treeDiscounts(2)/treeDiscounts(1);

end % function computation_forward_discount