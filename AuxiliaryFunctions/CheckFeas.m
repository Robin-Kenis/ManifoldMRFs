function B = CheckFeas(X, tol)
    % A simple metric to see how far a measure is from being a primal
    % feasible point
    B = zeros(size(X, 3), 1);
    for b = 1:size(X, 3)
        D = eig(X(:, :, b));
        B(b) = all(D >= -tol);
    end
end