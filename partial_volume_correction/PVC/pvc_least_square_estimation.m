function [ Resultats ] = pvc_least_square_estimation(intensite, X)
%
% Function that uses the least square estimation formula to solve a linear
% system 
%
% Resultats = pvc_least_square_estimation(intensite,X)
%
% Inputs:
%       intensite: Nx1 matrix containing the value of the mixed voxels. N
%           is the number of voxels
%       X: NxM mtrix containing the proportion of each tissue for each
%           voxel. N is the number of voxels and M the number of different
%           tissues
% Outputs:
%      Resultats: 2xM matrix containing in its first line the number of
%           the tissue and in its second line the value of the unmixed
%           tissue
%
% Function is written by Camille Van Assel Université de Montreal (november
% 2016)


Resultats = [1,2,3;0 0 0];

for iter=3:-1:1
    if sum(X(:,iter))==0
        X(:,iter)=[];
        Resultats(:,iter)=[];
    end
end

Resultats(2,:) = inv(X'*X)*X'*intensite;

end


