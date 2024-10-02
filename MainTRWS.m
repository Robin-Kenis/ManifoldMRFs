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

% Choose number of labels
NLabels = 150;

% Load problem instance
[Experiment, ~, ~, Image] = Experiments.Load(1);  % Load experiment 1, 2 or 3

% Define edge set of N1xN2-grid graph
E = EdgesGridGraph(Experiment.N1, Experiment.N2);

% Sample the Manifold
all_samples = 8*NLabels;
V = [];
num_samples = 0;
while num_samples ~= NLabels
    [V,Tri] = SpiralSampleSphere(all_samples, false);
    V = V(all(V > 0, 2), :);
    num_samples = size(V,  1);
    all_samples = all_samples + NLabels - num_samples;
end

% Construct the dataterms, regularization terms and adjacency matrix
NEdges = 2*Experiment.N1*Experiment.N2-Experiment.N1-Experiment.N2;
NVert = Experiment.N1*Experiment.N2;

CData   = reshape(permute(Image.ChromaticityNoisy, [2, 1, 3]), NVert, 3);
FV  = ones(size(V, 1), 1)*vecnorm(CData, 2, 2).' - V*CData.';
%FV  = FV(:);
%FV = reshape(FV, NLabels, NVert);

FE = Experiment.kconst*abs(acos(V*V.'));

adjacency = sparse(E(:, 1), E(:, 2), ones(NEdges, 1), NVert, NVert);

% Don't keep extra copy of large matrices
clear E CData num_samples;

%% Use TRW-S Solver to solve MRF

options.maxIter = 10000;
options.funcEps = 1e-4;
options.maxTime = Experiment.max_time;
options.verbosity = 2;

tic;
[labels, potential, lower_bound] = mrfMinimizeMex(FV, adjacency, FE, options);
time = toc;

%% Extract the solutions

% Extract solution marginals (Tmu_i, Tmu_ij)
vx          = V(labels, :).';

%% Check solution quality

optimal_lower_bound = lower_bound;
obtained_value      = potential;
gap                 = obtained_value - optimal_lower_bound;
fprintf('primal-dual gap %d \n', gap);
fprintf('primal %d \n', obtained_value);
fprintf('dual %d \n', optimal_lower_bound);


%% Figures

opt_chrom           = permute(reshape(vx.', Experiment.N1, Experiment.N2, 3), [2, 1, 3]);

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
