classdef S2_
    properties (Constant)
        NPVar               =   3;
        NDVar               =   3;
        Constraints         =   {[0,   0,   0,  -1;
                                 2,   0,   0,   1;
                                 0,   2,   0,   1;
                                 0,   0,   2,   1;];}
        ProjectionIndices   =  [ 0,  -3,   2;
                                 3,   0,  -1;
                                -2,   1,   0;];
    end
end