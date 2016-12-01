function VP  = pvc_cube_fraction( Volume, size_pixel,number_of_surf)

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
