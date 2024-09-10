function Basis = BorderBasisMatrices(Basis, A)
    % Takes in basis matrices and reduces them according to polynomial 
    % constraints in the matrix A
    for i = size(A, 2):-1:1
        Acol    = A(:, i);
        AMI     = find(Acol,1,'last');
        
        ANZ     = find(Acol);
        ANZ     = ANZ(1:end-1);
        
        for CMI = flip(ANZ.')
            Basis{CMI} = Basis{CMI} - A(CMI, i)*Basis{AMI}/A(AMI, i);
        end
        if AMI <= size(Basis{1})
            for l = 1:length(Basis)
                Basis{l}(AMI, :) = [];
                Basis{l}(:, AMI) = [];
            end
        end
        Basis(AMI)  = [];
        
    end
    
end