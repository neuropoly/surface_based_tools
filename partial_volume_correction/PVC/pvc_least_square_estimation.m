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
% 
% Function is written by Camille Van Assel Univertsité polytechnique de
% Montreal (november 2016). 
% 
% The MIT License (MIT)
% Copyright (c) 2016 Ecole Polytechnique, Université de Montréal
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


Resultats = [1,2,3;0 0 0];

for iter=3:-1:1
    if sum(X(:,iter))==0
        X(:,iter)=[];
        Resultats(:,iter)=[];
    end
end

Resultats(2,:) = inv(X'*X)*X'*intensite;

end


