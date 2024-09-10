function M = ConstructMM(mu, Basis)
    % Takes a basis for a LMI (in the form of a cell of matrices) and
    % a vector/matrix containing the weights of the LMI(s) as entries/rows
    % and returns the constructed LMI matrix/matrices
    M = zeros([size(Basis{1}), size(mu, 2)]);
    for i = 1:size(mu, 2)
        for j = 1:size(mu, 1)
            M(:, :, i) = M(:, :, i) + mu(j, i)*Basis{j};
        end
    end
end