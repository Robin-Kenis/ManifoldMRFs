classdef SO3_
    properties (Constant)
        % 9 primal variables (= dimension of the embedding space)
        NPVar               =   9;

        % 3 dual variables (= dimension of the tangent space)
        NDVar               =   3;

        % Cell containing encoded polynomial constraints p_i(x) = 0
        % Row format: [x_11, x_12, x_13, x_21, x_22, x_23, x_31, x_32, x_33, coefficient]
        % Every row encodes a monomial, with entries encoding the power of
        % the variable x_ij appearing in the monomial, and the coefficient
        % given as the final entry. The sum of the monomials represented by
        % the rows should equal zero on the given manifold.
        %
        % Example:      the row
        %          [1, 0, 0,   2, 0, 0,   3, 0, 0,     4]
        % encodes the monomial
        %       (1*x_11)*(2*x_12)*(3*x13)*(4*1) = 4*x_11*x_12^2*x_13^3,
        % and the matrix
        %          [0, 0, 0,   0, 0, 0,   0, 0, 0,    -1]
        %          [1, 0, 0,   2, 0, 0,   3, 0, 0,     4]
        % encodes the polynomial constraint
        %       -1 + 4*x_11*x_12^2*x_13^3 = 0.
        Constraints         =  {[0, 0, 0, 0, 0, 0, 0, 0, 0, -1;     % Constraint 1
                                 2, 0, 0, 0, 0, 0, 0, 0, 0,  1;     %   -1 + x_11^2 + x12^2 + x13^2 = 0
                                 0, 2, 0, 0, 0, 0, 0, 0, 0,  1;
                                 0, 0, 2, 0, 0, 0, 0, 0, 0,  1];
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, -1;     % Constraint 2
                                 0, 0, 0, 2, 0, 0, 0, 0, 0,  1;     %   -1 + x_21^2 + x22^2 + x23^2 = 0
                                 0, 0, 0, 0, 2, 0, 0, 0, 0,  1;
                                 0, 0, 0, 0, 0, 2, 0, 0, 0,  1];
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, -1;     % Constraint 3
                                 0, 0, 0, 0, 0, 0, 2, 0, 0,  1;     %   -1 + x_31^2 + x32^2 + x33^2 = 0
                                 0, 0, 0, 0, 0, 0, 0, 2, 0,  1;
                                 0, 0, 0, 0, 0, 0, 0, 0, 2,  1];
                                [0, 1, 0, 0, 0, 1, 0, 0, 0,  1;     % Constraint 4
                                 0, 0, 1, 0, 1, 0, 0, 0, 0,  -1;    %   x12*x23 - x13*x22 - x_31 = 0
                                 0, 0, 0, 0, 0, 0, 1, 0, 0,  -1];
                                [1, 0, 0, 0, 0, 1, 0, 0, 0,  -1;    % Constraint 5
                                 0, 0, 1, 1, 0, 0, 0, 0, 0,  1;     %   -x_11*x_23 + x23*x21 - x32 = 0
                                 0, 0, 0, 0, 0, 0, 0, 1, 0,  -1];
                                [1, 0, 0, 0, 1, 0, 0, 0, 0,  1;     % Constraint 6
                                 0, 1, 0, 1, 0, 0, 0, 0, 0,  -1;    %   x_11*x_22 - x12*x21 - x33 = 0
                                 0, 0, 0, 0, 0, 0, 0, 0, 1,  -1];
                                [1, 0, 0, 1, 0, 0, 0, 0, 0,  1;     % Constraint 7
                                 0, 1, 0, 0, 1, 0, 0, 0, 0,  1;     %   x_11*x21 + x12*x22 + x13*x23 = 0
                                 0, 0, 1, 0, 0, 1, 0, 0, 0,  1];
                                [1, 0, 0, 0, 0, 0, 1, 0, 0,  1;     % Constraint 8
                                 0, 1, 0, 0, 0, 0, 0, 1, 0,  1;     %   x_11*x31 + x12*x32 + x13*x33 = 0
                                 0, 0, 1, 0, 0, 0, 0, 0, 1,  1];
                                [0, 0, 0, 1, 0, 0, 1, 0, 0,  1;     % Constraint 9
                                 0, 0, 0, 0, 1, 0, 0, 1, 0,  1;     %   x_21*x31 + x22*x32 + x23*x33 = 0
                                 0, 0, 0, 0, 0, 1, 0, 0, 1,  1];}

        % Matrix representing the projection onto the manifold.
        ProjectionIndices   =  [ 0,  -3,   2;
                                 3,   0,  -1;
                                -2,   1,   0;
                                 0,  -6,   5;
                                 6,   0,  -4;
                                -5,   4,   0;
                                 0,  -9,   8;
                                 9,   0,  -7;
                                -8,   7,   0;];
    end
end