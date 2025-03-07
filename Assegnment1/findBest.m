function bestM = findBest(M, error, threshold)
% Find the best M 
%
%INPUT
% M:             vector of time steps for CRR/simulations
% error:         error vector of the relative value M
% threshold:     value of the bid/ask basis point

for i = 1:length(M)
    if (error(i) < threshold)
        bestM = M(i);
        break;
    end
end

end %function findBest