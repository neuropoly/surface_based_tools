function [ S ] = pvc_partial_volume_estimation(Vol, fsdir, voxel, output_path, nb_surface, second_correction)
% 
% This function pvc_compute_unmixed_values will recover the unmixed values
% of voxels prone to Partial Volume Effect. First step of the process will
% be to compute the partial volume fraction using the function
% pvc_compute_volume_fraction. Then, it will calculate the unmixed values
% using a linear regression with the neigbors voxels.
%
% S = pvc_compute_unmixed_values(Volume, fsdir, voxel, output_path, second_correction)
%
% Inputs:
%       Volume: the input volume that is prone to partial volume effect
%       (needs to be a .mgz or .mgh file)
%       fsdir: the subject directory
%       voxel: Resolution of the input volume
%            
%       (Optionnal) 
%       output_path : Outpout path
%       nb_surface : number of expanded surfaces created. The more surfaces, 
%           the more precise the method is. But the creation of surfaces 
%           is also very time consuming  
%       second_correction: Defaut false. Second correction for the partial
%       volume estimation. Can only be used in case of isotropic volume
%       
%
% Outputs:
% The function create several files in the current directory:
%           - vp_WM.mgz : mask of the white matter
%           - vp_CTX.mgz : mask of the cortex
%           - vp_CSF.mgz :  mask of the CSF
%           - vp_r_square.mgz:mask of the r square coefficients 
%           - VPC.mat: VPC is a structure containing: - the number of pixel crossed 
%                                         by the surface
%                               - the mean R-squared coefficient 
%                               - the % of R-squared coeeficient up to 0.8    
%                               - the % of R-squared coefficient up to 0.9
%                                       
% Function is written by Camille Van Assel Universit� polytechnique de
% Montreal (november 2016) - camille260395@gmail.com - gabriel.mangeat@gmail.com 

%% Check if there is optionnal arguments
if nargin < 4
    second_correction = false;
    output_path = pwd;
    nb_surface = 5;
end

if nargin ==4
    if isa(output_path,'logical')
        second_correction = true;
        output_path = pwd;
    elseif isa(output_path,'char')
        second_correction = false;
    elseif isa(output_path,'float')
        output_path = pwd;
        second_correction = false;
    else
        error('Error. \n Input must be a char (ouput path) or a logical (second correction) , not a %s.',class(output_path));
    end
end

if nargin == 5
    if isa(output_path,'float')
        second_correction = nb_surface;
        nb_surface = output_path;
    elseif isa(output_path,'char')
        if isa(nb_surface,'float')
            second_correction = false;
        elseif isa(nb_surface,'logical')
            second_correction = nb_surface;
            nb_surface = 5;
        end
    else
        error('Error. \n Input type is false');
    end 
end


%% First step: calculation of the fraction of each tissue in the pixels 
tic;

[~,~,ext] = fileparts(Vol);
if strcmp(ext ,'.nii') || strcmp(ext ,'.nii.gz')
    Mri_nii = load_nii(Vol);
    mri = Mri_nii.vol;
elseif strcmp(ext ,'.mgz') || strcmp(ext ,'.mgh')
    [mri,N] = load_mgh(Vol);
else
    error('Wrong input volume type. The volume has to be a .nii or .mgz file...');
end

%Creation of the expanded surfaces and computation of the partial
%volume fractions 
fsdir =[fsdir '/surf']; 

pvc_create_expanded_surfaces('lh.pial', fsdir, 5, voxel);
pvc_create_expanded_surfaces('lh.white',fsdir, 5, voxel);

VP_pial_lh = pvc_compute_volume_fraction('lh.pial', fsdir ,nb_surface, Vol,[output_path '/VP_pial_lh'],second_correction);
VP_white_lh = pvc_compute_volume_fraction('lh.white', fsdir ,nb_surface, Vol,[output_path '/VP_white_lh'],second_correction);

pvc_create_expanded_surfaces('rh.pial', fsdir, 5, voxel);
pvc_create_expanded_surfaces('rh.white',fsdir, 5, voxel);  

VP_pial_rh = pvc_compute_volume_fraction('rh.pial', fsdir ,nb_surface, Vol,[output_path '/VP_pial_rh'],'pial','rh',second_correction);
VP_white_rh = pvc_compute_volume_fraction('rh.white', fsdir ,nb_surface, Vol,[output_path '/VP_white_rh'],'white','rh',second_correction);  


%%  Second step: computation of the pixel value for each tissue

hemi = ['lh','rh'];
for iter=1:2
    pial_vol = ['VP_pial_' hemi(iter*2-1:iter*2)];
    WM_vol = ['VP_white_' hemi(iter*2-1:iter*2)];

    PIAL = eval(pial_vol);
    WM = eval(WM_vol);

    % Classe is a table that contains all the volume fraction 
    Classe = zeros(size(mri,1),size(mri,2),size(mri,3),3);

    Classe(:,:,:,2) = PIAL(:,:,:,1)-WM(:,:,:,1);
    Classe(:,:,:,1) = WM(:,:,:,1);
    Classe(:,:,:,3) = 1-PIAL(:,:,:,1);
    
    % Volume is a multilabel file containing the segmentation of the volume
    % Label 1 for the voxel prone to PVE, 2 for the voxel in the WM and
    % between 0 and 1 for the voxel prone to PVE
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

    for j=1:size(Resultats,1)
        img_1(Resultats(j,1),Resultats(j,2),Resultats(j,3)) = Moyenne(j,2);
        img_2(Resultats(j,1),Resultats(j,2),Resultats(j,3)) = Moyenne(j,1);
        img_3(Resultats(j,1),Resultats(j,2),Resultats(j,3)) = Moyenne(j,3);
    end
    
    if strcmp(ext ,'.nii') || strcmp(ext ,'.nii.gz')
        Mri_nii.vol=img_1;
        save_nifti(Mri_nii, [output_path '/vp_CTX_' side(iter*2-1:iter*2) '.nii.gz']);
        Mri_nii.vol=img_2;
        save_nifti(Mri_nii, [output_path '/vp_WM_' side(iter*2-1:iter*2) '.nii.gz']);
        Mri_nii.vol=img_3;
        save_nifti(Mri_nii, [output_path '/vp_CSF_' side(iter*2-1:iter*2) '.nii.gz']);
        Mri_nii.vol = R_square;
        save_nifti(Mri_nii, [output_path '/vp_r_square_' side(iter*2-1:iter*2) '.nii.gz']);
    else strcmp(ext,'.mgz') || strcmp(ext ,'.mgh') 
        save_mgh(img_1, [output_path '/vp_CTX_' side(iter*2-1:iter*2) '.mgz'],N);
        save_mgh(img_2, [output_path '/vp_WM_' side(iter*2-1:iter*2) '.mgz'],N);
        save_mgh(img_3, [output_path '/vp_CSF_' side(iter*2-1:iter*2) '.mgz'],N);  
        save_mgh(R_square, [output_path '/vp_r_square_' side(iter*2-1:iter*2) '.mgz'],N);
    end


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
    filename = ['VPC_' iter*2-1:iter*2 '.mat'];
    save(filename,'S');
end
toc;

end