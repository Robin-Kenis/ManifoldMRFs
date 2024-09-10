classdef Experiments < handle
    properties (Constant)
        % Imaging experiments, manifold is 'S2'
        Experiment = {struct(ImageFilename = 'Rooftop.jpg',     MaskFilename = 'Mask2.png',     D = 2,  N1 = 300,   N2 = 200,   max_time = 1000,    kconst = 3e-1,  sigma = 7.5e1);
                      struct(ImageFilename = 'Flowers.png',     MaskFilename = [],              D = 2,  N1 = 300,   N2 = 200,   max_time = 1000,    kconst = 3e-1,  sigma = 7.5e1);
                      struct(ImageFilename = 'Butterfly.jpg',   MaskFilename = [],              D = 2,  N1 = 300,   N2 = 200,   max_time = 1000,    kconst = 3e-1,  sigma = 7.5e1);}
    end
    methods (Static)
        function [Experiment, FV, FE, Image] = Load(index)            
            % Load image and extract intensity and chromaticity
            Experiment          = Experiments.Experiment{index};

            D                   = Experiment.D;
            N1                  = Experiment.N1;
            N2                  = Experiment.N2;
            kconst              = Experiment.kconst;
            sigma               = Experiment.sigma;

            DVecV               = nchoosek(2*ceil(D/2) + 3, 3);

            filename            = Experiment.ImageFilename;
            image_orig          = imread(filename);
            image_orig          = double(imresize(image_orig, [N2, N1]));
            intensity_orig      = vecnorm(double(image_orig), 2, 3);
            chromaticity_orig   = double(image_orig)./vecnorm(double(image_orig), 2, 3);
            
            % Add mask (for inpainting) and noise to image, and extract intensity and chromaticity
            maskname                = Experiment.MaskFilename;
            if ~isempty(maskname)
                mask                = permute(imread(maskname), [2, 1, 3]);
                if length(size(mask)) > 2
                    mask = rgb2gray(mask);
                end
                mask                = mask >= 100;
                mask                = double(imresize(mask, [N2, N1]));
            end
            image                   = image_orig + sigma*randn(size(image_orig)).*intensity_orig/256;
            if ~isempty(maskname)                   % This is statement should ideally have been 
                                                    % after the projection onto RGB, but for the 
                                                    % experiments in the paper it was accidentally
                                                    % put after (which is importantly an equally valid
                                                    % experiment)
            image(image < 0)        = 0;            % Project image back onto RGB-space
            image(image > 255)      = 255;          % Project image back onto RGB-space
                image               = image.*mask;
            end
            intensity               = vecnorm(double(image), 2, 3);
            chromaticity            = double(image)./vecnorm(double(image), 2, 3);
            
            % Setup cost functions
            if ~isempty(maskname)
                bound = 10;                         % Simulate that the mask is unknown;
                mask_est = intensity.';             % pixels with very low intensity are part of the estimated mask
                mask_est = (mask_est(:) < bound);   % Compute indices that are part of the mask
            end
            
            % compute matrix if coefficient vectors for vertex cost (data terms)
            CHROMdata           = reshape(permute(chromaticity, [2, 1, 3]), N2*N1, 3);    
            FV                  = sparse([1/2*ones(1, N1*N2); -CHROMdata.'; 1/2*ones(1, N1*N2); zeros(2, N1*N2); 1/2*ones(1, N1*N2); zeros(1, N1*N2); 1/2*ones(1, N1*N2)]);
            FV                  = [FV; sparse([], [], [], DVecV-10, N1*N2)];
            if ~isempty(maskname)
                FV(:, mask_est)     = 0;
            end
            
            % The regularization cost is the same (geodesic distance cost) for each
            % edge
            FE          = kron(ones(1, 2*N1*N2-N1-N2), [kconst; zeros((2*3+1)*DVecV - 1, 1)]);

            Image.ImageClean        = image_orig;
            Image.IntensityClean    = intensity_orig;
            Image.ChromaticityClean = chromaticity_orig;

            Image.ImageNoisy        = image;
            Image.IntensityNoisy    = intensity;
            Image.ChromaticityNoisy = chromaticity;

            Image.DenoisedIntensity = intensity_orig;
        end
    end
end