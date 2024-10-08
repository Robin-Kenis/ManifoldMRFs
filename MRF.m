classdef MRF < handle

    properties

        %% Manifold
        Manifold;           % Should be either "S1", "S2" or "SO3"
        ProjectionMatrix;   % Matrix projecting a polynomial vector v onto the tangent vector space of Manifold
        NPVar;              % Number of primal variables (dimension of the embedding space)
        NDVar;              % Number of dual variables (equal to or less than NPVar, depending on the projection matrix)

        %% Vertices
        NV;                 % Number of vertices of the MRF graph
        DV;                 % Degree of measure vector for every vertex
        FV;                 % Matrix of coefficient vectors encoding the data terms
        DVecV;              % Dimension of the measure vector for every vertex
        DMatV;              % Dimension of the measure matrix for every vertex
        BMatV;              % Basis matrices for PSD constraint for every vertex
        BVecV;              % Vectorized basis matrices for PSD constraint for every vertex 
        CstMV;              % Constraint matrix for measure vector for every vertex

        %% Edges
        E;                  % List of all the edges (as vertex pairs, enumerated 1 to NV)
        NE;                 % Number of edges of the MRF graph
        DE;                 % Degree of measure vector for every edge
        FE;                 % Matrix of coefficient vectors encoding the regularization terms
        DVecE;              % Dimension of the measure vector for every edge
        DMatE;              % Dimension of the measure matrix for every edge
        BMatE;              % Basis matrices for PSD constraint for every edge
        BVecE;              % Vectorized basis matrices for PSD constraint for every edge 
        CstME;              % Constraint matrix for measure vector for every vertex

    end
    methods
        function obj = MRF(Manifold, D, E)
            addpath("Manifolds\")
            addpath("AuxiliaryFunctions\")
            if isstring(Manifold)
                obj.Manifold    = eval(strcat(Manifold, "_"));
            else 
                obj.Manifold    = Manifold;
            end
            
            obj.NV          = max(E, [], 'all');
            obj.DV          = D;

            obj.E           = E;
            obj.NE          = size(E, 1);
            obj.DE          = D;
            
            obj.SetupMRF();
        end

        % Some functionality to setup the MRF structure, should be called
        % to properly initialize the MRF
        function obj = SetupMRF(obj)
            try
                obj.NPVar               = obj.Manifold.NPVar;
                obj.NDVar               = obj.Manifold.NDVar;

                obj.DMatV   = nchoosek(ceil(obj.DV/2) + obj.NPVar, obj.NPVar);
                obj.DMatE   = 2*obj.DMatV;
                obj.DVecV   = nchoosek(2*ceil(obj.DV/2) + obj.NPVar, obj.NPVar);
                obj.DVecE   = 3*obj.DVecV;

                ProjMatCell             = ConstructTangentProjectionMatrix(obj.NPVar, obj.NDVar, obj.DV, obj.Manifold.ProjectionIndices);
                obj.ProjectionMatrix    = sparse([], [], [], (obj.DVecV-1), (1+2*obj.NDVar)*obj.DVecV);
                for dv = 1:obj.NDVar
                    obj.ProjectionMatrix(1:obj.DVecV-1, 1+(1+2*(dv-1))*obj.DVecV:(2*dv)*obj.DVecV) = ProjMatCell{dv}(2:end, :);
                end

                obj.CstMV = [];
                for i = 1:length(obj.Manifold.Constraints)
                    obj.CstMV = [obj.CstMV, OrderedConstraintMatrix(obj.Manifold.Constraints{i}, obj.DV, obj.NPVar)];
                end

                obj.CstME = [];
            catch
                warning(strcat("The manifold is not known, or not properly configured. The setup could not be completed."));
                return;
            end

            % Setup constraint matrix
            obj.BMatV     = cell(obj.DVecV, 1);
            for p = 1:obj.DVecV
                deg             = IndexToDegree(p, obj.NPVar);
                obj.BMatV{p}    = construct_monomial_matrices(deg, obj.DV, obj.NPVar);
            end

            obj.EdgeBasisForVertexBasis();
            obj = VectorizeBasisMatrices(obj);
            
            ProjMatCell    = ConstructTangentProjectionMatrix(obj.NPVar, obj.NDVar, obj.DV, obj.Manifold.ProjectionIndices);
            obj.ProjectionMatrix     = sparse([], [], [], (obj.DVecV-1), (1+2*obj.NDVar)*obj.DVecV);
            for dv = 1:obj.NDVar
                obj.ProjectionMatrix(1:obj.DVecV-1, 1+(1+2*(dv-1))*obj.DVecV:(2*dv)*obj.DVecV) = ProjMatCell{dv}(2:end, :);
            end
        end

        function obj = ReduceToQuotientRing(obj)
            MatRedV     = ConstraintsToReductionMatrix(obj.CstMV);
            MatRedE     = kron(eye(2*obj.NDVar+1), obj.CstMV);
            
            if isempty(obj.FV)
                disp("The data terms are not yet specified, you should instantiate FV first.");
            else
                obj.FV = ReduceToBB(obj.FV, MatRedV);
            end
            if isempty(obj.FE)
                disp("The regularization terms are not yet specified, you should instantiate FE first.");
            else
                obj.FE = ReduceToBB(obj.FE, MatRedE);
            end

            obj.BMatV   = BorderBasisMatrices(obj.BMatV, MatRedV);

            obj.DVecV   = length(obj.BMatV);
            obj.DMatV   = size(obj.BMatV{1}, 1);
            obj.DMatE   = 2*obj.DMatV;
            obj.DVecE   = 3*obj.DVecV;

            obj.EdgeBasisForVertexBasis();
            obj = VectorizeBasisMatrices(obj);

            ProjMatCell    = ConstructTangentProjectionMatrix(obj.NPVar, obj.NDVar, obj.DV, obj.Manifold.ProjectionIndices);
            for j = 1:obj.NDVar
                ProjMatCell{j} = ReduceToBB(ProjMatCell{j}, MatRedV).';
                ProjMatCell{j} = ReduceVandermonde(ProjMatCell{j}, MatRedV);
            end
            obj.ProjectionMatrix     = sparse([], [], [], (obj.DVecV-1), (1+2*obj.NDVar)*obj.DVecV);
            for dv = 1:obj.NDVar
                obj.ProjectionMatrix(1:obj.DVecV-1, 1+(1+2*(dv-1))*obj.DVecV:(2*dv)*obj.DVecV) = ProjMatCell{dv}(2:end, :);
            end
        end

        function obj = EdgeBasisForVertexBasis(obj)
            obj.BMatE = cell(3*length(obj.BMatV), 1);
            for i = 1:2
                for j = i:2
                    for k = 1:length(obj.BMatV)
                        T = sparse([], [], [], size(obj.BMatV{1}, 1)*2, size(obj.BMatV{1}, 2)*2);
                        T(1+(i-1)*size(obj.BMatV{1}, 1):i*size(obj.BMatV{1}, 1), 1+(j-1)*size(obj.BMatV{1}, 2):j*size(obj.BMatV{1}, 2)) = obj.BMatV{k};
                        if i ~=j
                            T = T + T.';
                        end
                        obj.BMatE{k+((i-1)*2+(j-1)-(i-1)*i/2)*length(obj.BMatV)} = T;
                    end
                end
            end
        end

        function obj = DataTerms(obj, FV)
            obj.FV = FV;
        end

        function obj = RegularizationTerms(obj, FE)
            obj.FE = FE;
        end

        function obj = CostFunctions(obj, FV, FE)
            obj.DataTerms(FV);
            obj.RegularizationTerms(FE);
        end

        function obj = VectorizeBasisMatrices(obj)
            %% Vectorized basis matrices
            LI      = [zeros(1, obj.DMatV); sqrt(2)*speye(obj.DMatV)];
            RI      = sparse([], [], [], 1, obj.DMatV);
            SI      = [];
            for i = 1:obj.DMatV
                LI          = LI(2:end, :);
                LI(1, i)    = 1;
                RI(i)       = 1;
                SI          = [SI; kron(LI, RI)];
                RI(i)       = 0;
            end
            
            LIJ     = [zeros(1, obj.DMatE); sqrt(2)*speye(obj.DMatE)];
            RIJ     = sparse([], [], [], 1, obj.DMatE);
            SIJ     = [];
            for i = 1:obj.DMatE
                LIJ         = LIJ(2:end, :);
                LIJ(1, i)   = 1;
                RIJ(i)      = 1;
                SIJ         = [SIJ; kron(LIJ, RIJ)];
                RIJ(i)      = 0;
            end
            
            % CBasis = Chebyshev_basis_matrices(Bmats_ij);
            
            BVecV_      = cellfun( @(X) SI*X(:),  obj.BMatV,  'UniformOutput', false);
            BVecE_      = cellfun( @(X) SIJ*X(:), obj.BMatE, 'UniformOutput', false);
            
            obj.BVecV   = [BVecV_{:}];
            
            BVecE_      = [BVecE_{:}];
            obj.BVecE   = sparse([], [], [], 2*obj.DMatV*(2*obj.DMatV+1)/2, (1+2*obj.NDVar)*obj.DVecV);
            for dv = 1:obj.NDVar
                obj.BVecE(1+(dv-1)*size(BVecE_, 1):dv*size(BVecE_, 1), 1:obj.DVecV) = BVecE_(:, 1:obj.DVecV);
                obj.BVecE(1+(dv-1)*size(BVecE_, 1):dv*size(BVecE_, 1), 1+(1+2*(dv-1))*obj.DVecV:2*dv*obj.DVecV) = BVecE_(:, 1+obj.DVecV:2*obj.DVecV);
                obj.BVecE(1+(dv-1)*size(BVecE_, 1):dv*size(BVecE_, 1), 1+2*dv*obj.DVecV:(1+2*dv)*obj.DVecV) = BVecE_(:, 1+2*obj.DVecV:3*obj.DVecV);
            end
            clear BVecE_;
        end

    end
end