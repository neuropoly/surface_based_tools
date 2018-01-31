function [ Moyenne,  R_square, X_i] = pvc_compute_unmixed_values( Resultats, mri, Volume, Classe)
%
% Function that computes the unmixed value of every crossed voxel. For each
% voxel intersected by the surface, the fraction of tissues in the
% neighboring voxels are taken into account to solve the linear system with
% a least square estimation. 
% Then, the R-square coeeficient is calculted for each linear regression. 
%
% [Moyenne, R_square, X_i] = pvc_compute_unmixed_values(Resultats, mri, 
% Volume, Classe);
%
% Inputs:
%       Resultats: A Nx3 matrix containing the coordinates of the
%               intersected voxels. N is the number of interseted voxels.
%       mri: Matrix of the input volume
%       Volume: Matrix containing the label of the differenet tissue. For 
%               example 0 for tissue 1, 1 for a crossed voxel and 2 for 
%               tissue 2 
%       Classe: 4 dimensions matrix containing the proportion of each
%           tissue. For example, the matrix Classe(:,:,:,1) contains the 
%           volume fraction of tissue 1
%
% Outputs:
%       Moyenne: A NxM matrix containing the unmixed value of each tissue
%               for each intersected voxel. N is the number of intersected 
%               voxel and M the number of different tissue
%       R-squared: Matrix containing the R-squared coefficient of each
%               linear regression. This matrix has the same size as the
%               input volume.
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


nb_classe = size(Classe,4);

% For each pixel, we take the neighbors pixel (not crossed) and we apply
% the least square estimation
h = waitbar(0,'Please wait...');
steps = size(Resultats,1);
Moyenne = zeros(size(Resultats,1),nb_classe);

X_i = zeros(size(Resultats,1),nb_classe);
R_square =zeros(size(mri));
% For each pixel, we extract the values of the neighbors, and we apply the
% max-likelihood algorithm to these values
for iter=1:size(Resultats,1)
    X = [];
    intensite = [];
    size_c = 1;
    bool = 1;
    
    % The values in the cube must have a minimum number of non-crossed
    % pixel (here the minimun is 3 for each tissue)
    while bool
        carre = mri(max(1,Resultats(iter,1)-size_c):min(size(mri,1),Resultats(iter,1)+size_c),max(1,Resultats(iter,2)-size_c):min(size(mri,2),Resultats(iter,2)+size_c),...
            max(1,Resultats(iter,3)-size_c):min(size(mri,3),Resultats(iter,3)+size_c));
        classe_carre = Volume(max(1,Resultats(iter,1)-size_c):min(size(mri,1),Resultats(iter,1)+size_c),max(1,Resultats(iter,2)-size_c):min(size(mri,2),Resultats(iter,2)+size_c),...
            max(1,Resultats(iter,3)-size_c):min(size(mri,3),Resultats(iter,3)+size_c));
        [l,c,p] = meshgrid(1:size(carre,1),1:size(carre,2),1:size(carre,3));
        
        A = find(classe_carre == 0);
        B = find(classe_carre == 2);
        bool=0;
        
        if size(A,1)<4 && size(B,1)<4
            bool =1;
            size_c = size_c + 1;
        end
        
        C = find(classe_carre ==1);
        if size(C)<3
            bool=1;
            size_c = size_c + 1;   
        end
    end
    
    % then we extract the value of the pixels and their proportion of each
    % tissue
    for i=1:size(carre,1)*size(carre,2)*size(carre,3)
        pos = [max(1,Resultats(iter,1)-size_c)+c(i)-1,max(1,Resultats(iter,2)-size_c)+l(i)-1,max(1,Resultats(iter,3)-size_c)+p(i)-1];
        X(end+1,:) = [Classe(pos(1),pos(2),pos(3),1),Classe(pos(1),pos(2),pos(3),2),Classe(pos(1),pos(2),pos(3),3)];
        intensite(end+1,:) = carre(c(i),l(i),p(i));
    end  
    
    % We compute the Max-likelihood algorithm 
    M = pvc_least_square_estimation(intensite,X);
    for i=1:size(M,2)
        Moyenne(iter,M(1,i)) = M(2,i);
    end
   
   X_i(iter,:) = [Classe(Resultats(iter,1),Resultats(iter,2),Resultats(iter,3),1), Classe(Resultats(iter,1),Resultats(iter,2),Resultats(iter,3),2),...
       Classe(Resultats(iter,1),Resultats(iter,2),Resultats(iter,3),3)];
  
   % Estimation of the error
   valeur_estimee = X*Moyenne(iter,:)';
   moyenne = mean(intensite)*ones(size(valeur_estimee));
   
   % Estimation of the r-square value
   R_square(Resultats(iter,1),Resultats(iter,2),Resultats(iter,3)) = max(0,1 - sum((intensite - valeur_estimee).^2)/sum((intensite - moyenne).^2));
   waitbar(iter/ steps)
   
   if Classe(Resultats(iter,1),Resultats(iter,2),Resultats(iter,3),3) == 0
        Moyenne(iter,3)=0;
    end
   if Classe(Resultats(iter,1),Resultats(iter,2),Resultats(iter,3),1) == 0
        Moyenne(iter,1)=0;
   end
end   

close(h)

end

