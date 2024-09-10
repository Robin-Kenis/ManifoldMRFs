function idxlist = FindIndexList(dvec, D)
    % Finds row and column indices for positions of monomial with
    % multidegree dvec to appear in monomial matrix
    dmat = SplitIndex(dvec, D);
    idxlist = zeros(size(dmat, 1), 1);
    for i = 1:size(dmat, 1)
          idxlist(i, 1) = DegreeToIndex(dmat(i, :), ceil(D/2));
          idxlist(i, 2) = DegreeToIndex(dvec - dmat(i, :), ceil(D/2));
    end
end