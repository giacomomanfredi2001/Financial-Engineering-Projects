function yearly_indexes  = ComputingYearlyIndexes(bermudanYearFrac, dt)
% Computation of grid indexes corresponding to the nearest swaps year fractions
%
% INPUTS:
% bermudanYearFrac:         year fractions of the swaps
% dt:                       grid step
   
    % Computation of the time grid of the tree
    T = bermudanYearFrac(end-1);
    grid = (0 : dt : T)';

    % Initialization of the vector of indexes
    yearly_indexes = zeros(length(bermudanYearFrac) - 2, 1);

    % Index computation
    for i = 1 : length(yearly_indexes)
        index = find(grid > bermudanYearFrac(i), 1);

        difference1 = abs(grid(index) - bermudanYearFrac(i));
        difference2 = abs(grid(index - 1) - bermudanYearFrac(i));

        if difference2 >= difference1
            yearly_indexes(i) = index;
        else
            yearly_indexes(i) = index -1;
        end
    end
    
end % function ComputingYearlyIndexes