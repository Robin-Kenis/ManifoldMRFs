function F = ReduceVandermonde(F, A)
    for i = size(A, 2):-1:1
        Acol    = A(:, i);
        AMI     = find(Acol,1,'last');
        
        F(AMI, :) = [];
    end
end