function constraints = OrderedConstraintMatrix(eq, D, K)
    % Maps a polynomial equality constraint to a linear matrix constraint
    % encoding that constraint for the moment vectors
    Deq         = max(sum(eq(:, 1:K).'));
    Pr          = nchoosek(2*ceil(D/2)-Deq+K, K);
    P           = nchoosek(2*ceil(D/2)+K, K);
    constraints = [];
    for p = 1:Pr
        deg     = IndexToDegree(p, K);
        pol     = eq(:, 1:K) + deg;
        index   = zeros(size(pol, 1), 1);
        for j = 1:size(pol, 1)
            index(j, 1) = DegreeToIndex(pol(j, :), D);
        end
        constraints = [constraints, sparse(index, 1, eq(:, K+1), P, 1)];
    end
end