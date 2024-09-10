function prob = setup_mosek_geodesic_problem(mrf)
    
    %% This code assumes that the MRF formulation has been reduced (mrf.ReduceToQuotientRing())
    [r, res] = mosekopt('symbcon');
    symbcon = res.symbcon;
    
    %% Cost function
    % Cost in terms of linearly constrained variables
    prob.c          = [mrf.FV(:); mrf.FE(:)];
    
    % Cost in terms of SDP constrained variables
    prob.barc.subj = [];
    prob.barc.subk = [];
    prob.barc.subl = [];
    prob.barc.val  = [];
    
    %% Equality constraints
    % Equality and inequality constraints -- linearly constrained variables
    NMconst_i       = mrf.DMatV*(mrf.DMatV+1)/2;
    NMconst_ij      = 2*mrf.DMatV*(2*mrf.DMatV+1)/2;
    
    asize           = [0, mrf.DVecV*mrf.NV + 2*(1+2*mrf.NDVar)*mrf.DVecV*mrf.NV];
    
    prob.a          = sparse([], [], [], asize(1), asize(2));
    
    add_A           = [kron(speye(mrf.NV), [1, sparse([], [], [], 1, mrf.DVecV-1)]), sparse([], [], [], mrf.NV, mrf.NE*(2*mrf.NDVar+1)*mrf.DVecV)];
    
    prob.a          = [prob.a; add_A];
    
    TV_constr       = [speye(mrf.DVecV), kron(ones(1, mrf.NDVar), [zeros(mrf.DVecV, mrf.DVecV), -speye(mrf.DVecV)])];
    add_A           = [sparse([], [], [], mrf.NE*mrf.DVecV, mrf.NV*mrf.DVecV), kron(speye(mrf.NE), TV_constr)];
    prob.a          = [prob.a; add_A];
    
    Sel             = [sparse([], [], [], mrf.DVecV-1, 1), speye(mrf.DVecV-1, mrf.DVecV-1)];
    OSS             = ones(mrf.NE, 1)*[1, -1];
    ISS             = (1:mrf.NE).'*[1, 1];
    SPS             = sparse(ISS(:), mrf.E(:), OSS(:), mrf.NE, mrf.NV);
    OPSPS1          = kron(SPS, Sel);
    OPSPS2          = kron(speye(mrf.NE), mrf.ProjectionMatrix);
    OPSPS           = [OPSPS1, OPSPS2];
    
    prob.a          = [prob.a; OPSPS];
    
    prob.blc        = sparse(1:mrf.NV, ones(mrf.NV, 1), ones(mrf.NV, 1), size(prob.a, 1), 1);
    prob.buc        = sparse(1:mrf.NV, ones(mrf.NV, 1), ones(mrf.NV, 1), size(prob.a, 1), 1);
    
    %% SDP constraints
    prob.f          = [kron(speye(mrf.NV), mrf.BVecV), sparse([], [], [], mrf.NV*NMconst_i, (2*mrf.NDVar+1)*mrf.DVecV*mrf.NE); sparse([], [], [], mrf.NDVar*mrf.NE*NMconst_ij, mrf.NV*mrf.DVecV), kron(speye(mrf.NE), mrf.BVecE)];
    prob.g          = [zeros(mrf.NV*NMconst_i+mrf.NE*mrf.NDVar*NMconst_ij, 1)];      
    
    prob.accs       = [kron(ones(1, mrf.NV), [symbcon.MSK_DOMAIN_SVEC_PSD_CONE NMconst_i]), kron(ones(1, mrf.NDVar*mrf.NE), [symbcon.MSK_DOMAIN_SVEC_PSD_CONE NMconst_ij])];
 
end

