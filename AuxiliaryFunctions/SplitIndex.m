function dmat = SplitIndex(dvec, D)
    % Takes an index according to lexicographical ordering of monomials and
    % splits it into row and column index for position of corresponding 
    % monomial in moment matrix
    if length(dvec) ~= size(dvec, 2)
        dvec = dvec.';
    end

    Dhalf = ceil(D/2);
    Ldvec = length(dvec);
    for i = Ldvec:-1:1
        Dcur        = min(D, dvec(i));
        O{i}        = ones(Dcur+1, Ldvec);
        O{i}(:, i)  = (0:Dcur).';
    end

    dmat = kr(O);
    dmat = -dmat + dvec;
    dmat = dmat(all(dmat >= 0, 2), :);
    dmat = dmat(sum(dmat, 2) <= Dhalf, :);
    dmat = -dmat + dvec;
    dmat = dmat(sum(dmat, 2) <= Dhalf, :);
end