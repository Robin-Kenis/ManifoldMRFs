classdef S1_
    properties (Constant)
        NPVar               =   2;
        NDVar               =   1;
        Constraints         =   {[  0,   0,  -1;
                                   2,   0,   1;
                                   0,   2,   1;
                                   0,   0,   1;]};
        ProjectionIndices   =   [ -2,   1;];
    end
end