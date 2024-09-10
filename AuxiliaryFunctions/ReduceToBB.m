function F = ReduceToBB(F, A)
    for i = size(A, 2):-1:1
        Acol    = A(:, i);
        ANZ     = find(Acol);

        AMI     = ANZ(end);
        ANZ     = ANZ(1:end-1);
        
        for CMI = flip(ANZ.')
            F(CMI, :) = F(CMI, :) - A(CMI, i)*F(AMI, :)/A(AMI, i);
        end
        F(AMI, :) = [];
    end
end