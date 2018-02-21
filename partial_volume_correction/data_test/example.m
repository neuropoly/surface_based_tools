%%% example

clc
close all
clear all

%%% set paraneters

% fsdir: here refers to the actual fsdir/Subject_name
% the script will assume that the surfaces (xh.white, xh.pial)
% are in fsdir/surf/
fsdir = 'YOUR_PATH/partial_volume_correction/data_test/ctrl_test';

% hemisphere to consider ('lh' or 'rh')
hemi = 'lh';

% Volume that is subject to partial volume
Vol = [fsdir '/mri/T1.mgz'];

% Isotropic voxel resolution in mm
voxel = 1; % 1 mm isotroptic

% ouptut directory
output_path = [fsdir '/output'];
mkdir(output_path);

% number of surfaces that will be extended, adding more 
% surface will give a more precise estimation but will increase the
% computation time
nb_surface = 5; 

% second correction taking into account the normal of the surface
second_correction = false;





%%% Example 1:
% run the main function: pvc_partial_volume_estimation

compute_all = 0;
if compute_all
    extend = true;
    pvc_partial_volume_estimation(Vol, fsdir, hemi, voxel, output_path, nb_surface, extend, second_correction)
end




%%% Example 2:
% Run step by step (for the left hemisphere):

% Step 1, expending surfaces
% create #nb_surface expended surfaces in the positive and negative
% direction, for pial amd white surfaces; so 2* 2* #nb_surface
% surf created in [fsdir '/surf'], e.g. lh.pial0000, lh.pial0001,
% lh.pialn0001, etc.
% Long step, ~55min per surf, so ~20h on a single i7 processor.

compute_expended_surfaces = 0;
if compute_expended_surfaces
    pvc_create_expanded_surfaces('lh.pial', fsdir, nb_surface, voxel);
    pvc_create_expanded_surfaces('lh.white',fsdir, nb_surface, voxel);
end

% Step 2, compute partial volumes and r-square coefficients of fitting
% ~30 min to compute on a single i7 processor

compute_partial_volumes = 0;
if compute_partial_volumes
    extend = false;
    pvc_partial_volume_estimation(Vol, fsdir, hemi, voxel, output_path, nb_surface, extend, second_correction)
end



