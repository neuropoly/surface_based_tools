function [ S ] = pvc_partial_volume_estimation(Vol, fsdir, voxel, surf,  output_path, aseg, reg, second_correction)
% 
% This function pvc_compute_unmixed_values will recover the unmixed values
% of voxels prone to Partial Volume Effect. First step of the process will
% be to compute the partial volume fraction using the function
% pvc_compute_volume_fraction. Then, it will calculate the unmixed values
% using a linear regression with the neigbors voxels.
%
% S = pvc_compute_unmixed_values(Volume, fsdir, voxel, cortical_study, aseg, type, reg)
%
% Inputs:
%       Volume: the input volume that is prone to partial volume effect
%       (needs to be a .mgz or .mgh file)
%       fsdir: the subject directory
%       voxel: Resolution of the input volume
%       output_path : Outpout path
%       
%       (Optionnal) 
%       aseg: in the case of the cortical study, the aseg file needs to be
%               given in order to fill the GM and WM volume.
%       reg: In case the aseg file is not the same size/resolution of the
%               input volume, the reg file needs to be given.
%       second_correction: Defaut false. Second correction for the partial
%       volume estimation. Can only be used in case of isotropic volume
%       
%
% Outputs:
% The function create several files in the current directory:
%           - vp_WM.mgz : mask of the white matter
%           - vp_CTX.mgz : mask of the cortex
%           - vp_CSF.mgz :  mask of the CSF
%           - VPC.mat: VPC is a structure containing: - the number of pixel crossed 
%                                         by the surface
%                               - the mean R-squared coefficient 
%                               - the % of R-squared coeeficient up to 0.8    
%                               - the % of R-squared coefficient up to 0.9
%                                       
% Function is written by Camille Van Assel Université polytechnique de
% Montreal (november 2016)

if nargin < 5 || nargin == 6
    second_correction = false;
end

tic;
nb_surface = 5;
[mri,N] = load_mgh(Vol);


%% First step: calculation of the fraction of each tissue in the pixels 

%Creation of the expanded surfaces and computation of the partial
%volume fractions 
fsidr =[fsdir '/surf']; 
pvc_create_expanded_surfaces('lh.pial', fsdir, 5, voxel);
pvc_create_expanded_surfaces('lh.white',fsdir, 5, voxel);

VP_pial = pvc_compute_volume_fraction('lh.pial', fsdir ,nb_surface, Vol, 'VP_pial.mgz', cortical_study, 'pial', aseg,reg,'lh');
VP_white = pvc_compute_volume_fraction('lh.white', fsdir ,nb_surface, Vol, 'VP_white.mgz', cortical_study, 'white', aseg,reg,'lh');

pvc_create_expanded_surfaces('rh.pial', fsdir, 5, voxel);
pvc_create_expanded_surfaces('rh.white',fsdir, 5, voxel);  

VP_pial = pvc_compute_volume_fraction('rh.pial', fsdir ,nb_surface, Vol, size(mri), second_correction, cortical_type, 'pial', aseg,reg,'rh');
VP_white = pvc_compute_volume_fraction('rh.white', fsdir ,nb_surface, Vol, size(mri), second_correction, cortical_type, 'white', aseg,reg,'rh');  


%%  Second step: computation of the pixel value for each tissue

pial_vol = strcat(current_dir, '/VP_pial.mgz');
WM_vol = strcat(current_dir, '/VP_white.mgz');

PIAL = load_mgh(pial_vol);
WM = load_mgh(WM_vol);

Classe = zeros(size(mri,1),size(mri,2),size(mri,3),3);

Classe(:,:,:,2) = PIAL(:,:,:,1)-WM(:,:,:,1);
Classe(:,:,:,1) = WM(:,:,:,1);
Classe(:,:,:,3) = 1-PIAL(:,:,:,1);
Volume = Classe(:,:,:,2);

% Assign the label 2 to the pixel which are in the white matter but not in
% the cortex
Volume(Classe(:,:,:,3)==1)=2;   


% We determine which pixels are crossed
level = [0.01, 0.99];
nb = imquantize(Volume, level);
    
% We extract the coordinates of the crossed pixel 
[L,C,P] = meshgrid(1:size(Volume,2),1:size(Volume,1),1:size(Volume,3));
x = find(nb==2);
Resultats = [C(x), L(x), P(x)];

% For each pixel, we take the neighbors pixel (not crossed) and we apply
% the least square estimation

[Moyenne, R_square] = pvc_compute_unmixed_values( Resultats, mri, Volume, Classe);

%% Step 3: Creation of the new images

img_1 = mri;
img_2 =mri;
img_3 = zeros(size(mri));

img_1(Classe(:,:,:,2) == 0) = 0;
img_2(Classe(:,:,:,1) == 0) = 0;

for iter=1:size(Resultats,1)
    img_1(Resultats(iter,1),Resultats(iter,2),Resultats(iter,3)) = Moyenne(iter,2);
    img_2(Resultats(iter,1),Resultats(iter,2),Resultats(iter,3)) = Moyenne(iter,1);
    if cortical_study
        img_3(Resultats(iter,1),Resultats(iter,2),Resultats(iter,3)) = Moyenne(iter,3);
    end
end

save_mgh(img_1, [outpout_path '/vp_CTX.mgz'],N);
save_mgh(img_2, [output_path '/vp_WM.mgz'],N);
save_mgh(img_3, [output_path '/vp_CSF.mgz'],N);  
save_mgh(R_square, [output_path '/vp_r_square.mgz'],N);


%% Step 4 : Estimation of the robustness of the method
% Calculation of parameters

B = (Classe(:,:,:,2)~=1 & Classe(:,:,:,2)~=0);
nb_coupe = sum(sum(sum(B)));

C = (R_square~=0);
moyenne_R = sum(sum(sum(R_square)))/sum(sum(sum(C)));
D = R_square > 0.9;
R_sup09 = sum(sum(sum(D)));

E = R_square>0.8;
R_sup08 = sum(sum(sum(E)));

S = struct('nb_coupe',nb_coupe,'moyenne_R',moyenne_R,'R_sup09',R_sup09,'R_sup08',R_sup08);
filename = 'VPC.mat';
save(filename,'S');
toc;

end