function [VP] = pvc_compute_volume_fraction(Vol, surf_name, surf_path, number_of_surf, image, output_file, second_correction)
%
% Function that determines the volume fraction of each pixel crossed by the
% surface that delineates the different tissues.
%
% VP = pvc_compute_volume_fraction(surf_name,surf_path,number_of_surf, 
% image,output_file,second_correction)
%
% Inputs:
%       surf_name: name of the surface delineating the tissues
%       surf_path: path of the directory containing the surface
%       number_of_surf: number of surface expanded on one side of the
%               surface
%       image: path of the input Volume we want to find the fvolume
%               fraction
%       output_file: name of the output file
%
%   (Optionnal)
%       Second_correction: Default False
%
%  Outputs:
%       - VP_out.mgz: Partial volume fraction image, created in the current
%               directory.
%
% Function is written by Camille Van Assel, Ecole Polytechnique, Univertsite de
% Montreal (november 2016). 
% 
% The MIT License (MIT)
% Copyright (c) 2016 Ecole Polytechnique, Univertsite de
% Montreal
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

if nargin < 7
    second_correction = false;
end

%Compute the volume fractions

[~,~,ext] = fileparts(Vol);
if strcmp(ext ,'.nii') || strcmp(ext ,'.nii.gz')
    Mri_nii = load_nii(Vol);
    mri = Mri_nii.vol;
elseif strcmp(ext ,'.mgz') || strcmp(ext ,'.mgh')
    [mri,M] = load_mgh(Vol);
else
    error('Wrong input volume type. The volume has to be a .nii or .mgz file...');
end

size_image = size(mri); 

% Get the name of the expanded surfaces
launch_bash_profile = 'source ~/.bash_profile;';
Volume = zeros(size_image(1),size_image(2),size_image(3),number_of_surf);

for iter=0:2*number_of_surf
    if iter<10 && iter<number_of_surf+1
        surf = [surf_path '/' surf_name '00' num2str(iter)];
    elseif iter<number_of_surf+1 && iter>=10
        surf = [surf_path '/' surf_name '0' num2str(iter)];
    elseif iter>number_of_surf && iter<number_of_surf+10
        surf = [surf_path '/' surf_name 'n00' num2str(iter-number_of_surf)];
    elseif iter>number_of_surf && iter>number_of_surf+10
        surf = [surf_path '/' surf_name 'n0' num2str(iter-number_of_surf)];
    end
    [vertex_coords, faces] = read_surf(surf);
    faces = faces+1;
    
    %% Find the offset between Surface RAS and RAS space
    offset = zeros(3,1);
    cmd = [launch_bash_profile 'mris_info ' surf];
    [~,res] = system(cmd);
    a = strfind(res,'c_(ras)');
    pos = a(1) + 11;
    
    l = ',,)';
    for j=1:3
        Pos=strfind(res(pos:end),l(j));
        offset(j,1) = str2double(res(pos:pos+Pos(1)-2));
        pos=Pos(1)+pos+1;
    end
    
    vertex_coords = vertex_coords + [offset(1,1)*ones(size(vertex_coords,1),1) offset(2,1)*ones(size(vertex_coords,1),1) offset(3,1)*ones(size(vertex_coords,1),1)];
    
    
    %% Convert FreeSurfer RAS values into voxel indices
    Nvertices = size(vertex_coords,1);
    right_column = [ ones( Nvertices, 1 ); 0 ];
    SurfVertices = [ [vertex_coords; 0 0 0]  right_column ];
    
    cmd = [launch_bash_profile 'mri_info ' image ' --ras2vox'];
    [~, res] = system(cmd);

    fsRAS2Vox = zeros(4,4);
    [l,c] = meshgrid(1:4,1:4);
    
    split = strsplit(res);
    for i=0:15
        fsRAS2Vox(l(16-i),c(16-i)) = str2double(split{1,size(split,2)-i-1});
    end
    
    VoxVertices = fsRAS2Vox * SurfVertices';
    VoxVertices = VoxVertices(1:3,:)';
    
    VoxVertices = VoxVertices(:,[2 1 3]);
    VoxVertices = max(VoxVertices, 0);
    
    %% Create a binary volume of the crossed voxel
    FV=struct('faces',faces,'vertices',VoxVertices);
    if iter ==0
        FV0 = FV;
    end
    Volume(:,:,:,iter+1) = polygon2voxel(FV,size_image,'none',false);
   
    
end

gridX = [0:size_image(1)-1];
gridY = [0:size_image(2)-1];
gridZ = [0:size_image(3)-1];

Fraction = VOXELISE(gridX,gridY,gridZ,FV0);
Fraction = permute(Fraction,[2 1 3]);

% 0-based matrix to 1-based matrix
Fraction = cat(1,Fraction(2:size(Fraction,1),:,:),Fraction(1,:,:));
Fraction = cat(2,Fraction(:,2:size(Fraction,2),:),Fraction(:,1,:));
Fraction = cat(3,Fraction(:,:,2:size(Fraction,3)),Fraction(:,:,1));

%% determine the fraction of partial volume
if second_correction
    VP = cube_fraction(Volume, pixel, number_of_surf);
else
    A = (Volume(:,:,:,1)~=0);
    S = sum(Volume,4);
    VP = A - A./(S+1);
    for i=1:number_of_surf
        B = (Volume(:,:,:,i+1)~=0 & A);
        VP = VP-B./(S+1);
        A = B;
    end
end

% Determine the maximal number of surfaces that can cut a cube, depending
% on its normal vector
% N = (Vecteurs_normaux~=0);
% Norm = sqrt(Vecteurs_normaux(:,:,:,1).^2+Vecteurs_normaux(:,:,:,2).^2+Vecteurs_normaux(:,:,:,1).^2);
% PS = (Vecteurs_normaux(:,:,:,1)./Norm).*sign(Vecteurs_normaux(:,:,:,1))+...
%     (Vecteurs_normaux(:,:,:,2)./Norm).*sign(Vecteurs_normaux(:,:,:,2))+...
%     (Vecteurs_normaux(:,:,:,3)./Norm).*sign(Vecteurs_normaux(:,:,:,3));
% VP(A) = sqrt(3)*VP(A)./PS(A);



%% Draw the pixel inside the surface in the case of cortical study

x = (Fraction~=0);
VP(VP==0 & x==1)=1;

% Since we had to add 1 to every vertice coordinates, the coordinates of the image if 1-based instead of 0-based
% We have to  delete the last ligne in each dimension and add a ligne at the begining of the matrix
VP = cat(1,VP(size(VP,1),:,:),VP(1:size(VP,1)-1,:,:));
VP = cat(2,VP(:,size(VP,2),:),VP(:,1:size(VP,2)-1,:));
VP = cat(3,VP(:,:,size(VP,3)),VP(:,:,1:size(VP,3)-1));

%% Save the image 

if strcmp(ext ,'.nii') || strcmp(ext ,'.nii.gz')
     name = [output_file '.nii.gz'];
     Mri_nii.vol = VP;
     save_nifti(Mri_nii,name);
else strcmp(ext ,'.mgz') || strcmp(ext ,'.mgh')
    name = [output_file '.mgz'];
    save_mgh(VP, name , M);
end

end
