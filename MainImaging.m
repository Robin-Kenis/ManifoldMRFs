%% Geodesic total variation denoising of S2-valued signals
%
%
%
%

clear all;
close all;
clc;

addpath('AuxiliaryFunctions');
rng(12345);

%% Graph structure and problem parameters

% Load problem instance
[Experiment, FV, FE, Image] = Experiments.Load(1);  % Load experiment 1, 2 or 3

% Define edge set of N1xN2-grid graph
E = EdgesGridGraph(Experiment.N1, Experiment.N2);

% Construct S2-valued MRF object, load cost and reduce problem
% dimensionality
mrf = MRF("S2", Experiment.D, E);
mrf.CostFunctions(FV, FE);
mrf.ReduceToQuotientRing();

% Don't keep extra copy of large matrices
clear E FV FE;

%% Setup optimization problem

prob = SetupMosekGeodesic(mrf);

%% Solve optimization problem

% Remove some unnecessary

param.MSK_DPAR_OPTIMIZER_MAX_TIME       = Experiment.max_time;
param.MSK_IPAR_INTPNT_MAX_ITERATIONS    = 50;
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS      = 1e-10;
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS      = 1e-10;
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP    = 1e-20;

disp('solver started');
tic;
[r,res] = mosekopt('minimize info', prob, param);
time    = toc;

%% Extract the solutions

% Extract solution marginals (Tmu_i, Tmu_ij)
xx          = res.sol.itr.xx;
Tmu_i       = reshape(xx(1:mrf.NV*mrf.DVecV), mrf.DVecV, mrf.NV);
Tmu_ij      = reshape(xx(1+mrf.NV*mrf.DVecV:end), (2*mrf.NDVar+1)*mrf.DVecV, mrf.NE);

% Matrices
X_i         = ConstructMM(Tmu_i, mrf.BMatV);
X_ij        = zeros(mrf.DMatE, mrf.DMatE, mrf.NE, mrf.NDVar);
for dv = 1:mrf.NDVar
    X_ij(:, :, :, dv) = ConstructMM(Tmu_ij([1:mrf.DVecV, (1+2*(dv-1))*mrf.DVecV+1:(1+2*dv)*mrf.DVecV], :), mrf.BMatE);
end

% Rounded solution (simply project first moments onto sphere)
vx          = Tmu_i(2:4, :).'./vecnorm(Tmu_i(2:4, :).', 2, 2);

%% Check solution quality

PV          = res.sol.itr.pobjval;
DV          = res.sol.itr.dobjval;

PoV         = full(sum(mrf.FV.*Tmu_i, 1));
PoE         = full(sum(mrf.FE.*Tmu_ij, 1));

f_v         = full(mrf.FV(1:4, :));
F_r         = sum(f_v.*[ones(1, size(vx, 1)); vx.'], 1);
F_s         = Experiment.kconst*abs(acos(sum(vx(mrf.E(:, 1), :).*vx(mrf.E(:, 2), :), 2)));

F_rSDP      = sum(PoV);
F_sSDP      = sum(PoE);

optimal_lower_bound = res.sol.itr.dobjval;
obtained_value      = sum(F_r) + sum(F_s);
gap                 = obtained_value - optimal_lower_bound;
fprintf('primal-dual gap %d \n', gap);
fprintf('primal %d \n', obtained_value);
fprintf('dual %d \n', optimal_lower_bound);

%% Check rank-1 condition

% A quick check for optimality of the vertex measures
Rvert = PrimFeas(X_i);
Bvert = CheckFeas(X_i, 1e-5);

%% Figures

opt_chrom           = permute(reshape(vx, Experiment.N1, Experiment.N2, 3), [2, 1, 3]);

image_orig          = uint8(round(Image.ImageClean));
image               = uint8(round(Image.ImageNoisy));
den_chromaticity    = uint8(round(256.*opt_chrom));
part_den_image      = uint8(round(Image.DenoisedIntensity.*Image.ChromaticityNoisy));
denoised_image      = uint8(round(Image.DenoisedIntensity.*opt_chrom));
figure();
subplot(2, 2, 1);
hold on;
imshow(image_orig);
title('reference image');
hold off;
subplot(2, 2, 2);
hold on;
imshow(image);
title('noisy image');
hold off;
subplot(2, 2, 3);
hold on;
imshow(den_chromaticity);
title('denoised chromaticity');
hold off;
subplot(2, 2, 4);
hold on;
imshow(denoised_image);
title('denoised image');
hold off;

%%

% cur_date = date();
% save(strcat("INPAINTING_", filename, "_", cur_date, "_D", num2str(DP), "_K", num2str(kconst), "_S", num2str(sigma), "_", num2str(N1), "_", num2str(N2), ".mat"));
