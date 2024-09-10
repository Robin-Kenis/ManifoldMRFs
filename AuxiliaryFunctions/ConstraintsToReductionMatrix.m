function A = ConstraintsToReductionMatrix(A)
    % Takes a matrix encoding polynomial equality constraints for the
    % (standard) basis of moment vectors and returns a matrix to be used in
    % computing a reduced basis of moment vectors
    if isempty(A)
        return
    end
    A = A.';
    A = A(:, end:-1:1);
    A = rref(A);
    A = A(:, end:-1:1);
    A = A.';
    A = A(:, end:-1:1);
    A(:, find(all(A==0))) = []; 
    A = sparse(A);
end