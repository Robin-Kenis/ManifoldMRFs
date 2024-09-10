function mtx = construct_monomial_matrices(deg, D, K)
    % Construct the basis of matrices that make up the moment LMI
    % constraint
    N           = nchoosek(ceil(D/2)+K, K);
    idxlist     = FindIndexList(deg, D);
    mtx         = sparse(idxlist(:, 1), idxlist(:, 2), ones(size(idxlist, 1), 1), N, N);
end