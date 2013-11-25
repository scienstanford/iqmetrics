function [M1,MW,MW1] = perform_steerable_matching(M1,M,options)

% M1 = perform_steerable_matching(M1,M,options);

% spatial equalization
M1 = perform_histogram_equalization(M1,M);
 
% forward transforms
Jmin = 4;
options.nb_orientations = 4;
MW1 = perform_steerable_transform(M1, Jmin, options);
if iscell(M)
    MW = M;
else
    MW = perform_steerable_transform(M, Jmin, options);
end
% wavelet domain equalization
MW1 = perform_histogram_equalization(MW1,MW);
% backward transform
M1 = perform_steerable_transform(MW1, Jmin, options);

% spatial equalization
M1 = perform_histogram_equalization(M1,M);