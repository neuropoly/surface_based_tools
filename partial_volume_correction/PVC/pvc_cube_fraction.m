function VP  = pvc_cube_fraction( Volume, size_pixel,number_of_surf)
% Internal function
% Compute fractions of cubes via the number of surfaces crossed and the
% normal orientation
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

A = (Volume(:,:,:,1)~=0);
S = sum(Volume,4);
VP = A - A./(S+2);
for i=1:number_of_surf
    B = (Volume(:,:,:,i+1)~=0 & A);
    VP = VP-B./(S+2);
    A = B;
end

%% case 1: the surface normal is parallel to one of the cube edge
alpha= VP;
alpha_lim = floor((2*nb_surface+1)/(3.2*size_pixel));
alpha(S>alpha_lim)=0;
alpha= cube_alpha(alpha,size_pixel);
VP(S<=alpha_lim)=alpha(S<=alpha_lim);

%% case 2: The surface normal is parallel to the first diagonal
beta = VP;
beta_lim = ceil((sqrt(2)*pixel_size*2*nb_surface+1)/(3.2*size_pixel));
beta(S<=alpha_lim|S>beta_lim)=0;
beta = cube_beta(beta,size_pixel);
VP(S>alpha_lim&S<=beta_lim)=beta (S>alpha_lim&S<=beta_lim);

%% case 3: The surface normal is parallel to the big diagonal
gamma = VP;
gamma(S<=beta_lim)=0;
gamma = cube_gamma(gamma,size_pixel);
VP(S>beta_lim)=gamma(S>beta_lim); 
end

function h = cube_alpha(h,a)
h = a^2*h;
end

function h = cube_beta(h,a)
h(h < a*sqrt(2)/2 & h>0)= a*h(h < a*sqrt(2)/2 & h>0).^2;
h(h > a*sqrt(2)/2)=a^3-a*(a*sqrt(2)-h(h > a*sqrt(2)/2)).^2;
end

function h=cube_gamma(h,a)
h(h<a*sqrt(3)/3 & h>0) = sqrt(3)/2*h(h<a*sqrt(3)/3 & h>0).^3;
h(h>a*sqrt(3)/3 & h<2*a*sqrt(3)/3) = 1/6*a^3 + sqrt(3)/2*a^2*(h(h>a*sqrt(3)/3 & h<2*a*sqrt(3)/3)-a*sqrt(3)/3) + 3/2*a*(h(h>a*sqrt(3)/3 & h<2*a*sqrt(3)/3)-sqrt(3)/3).^2 -sqrt(3)*(h(h>a*sqrt(3)/3 & h<2*a*sqrt(3)/3)-a*sqrt(3)/3).^3;
h(h>2*a*sqrt(3)/3) = a^3-sqrt(3)/2*(a*sqrt(3)-h(h>2*a*sqrt(3)/3)).^3;
end
