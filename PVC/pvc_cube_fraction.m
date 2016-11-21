function VP  = pvc_cube_fraction( Volume, size_pixel,number_of_surf)

A = (Volume(:,:,:,1)~=0);
S = sum(Volume,4);
VP = A - A./(S+2);
for i=1:number_of_surf
    B = (Volume(:,:,:,i+1)~=0 & A);
    VP = VP-B./(S+2);
    A = B;
end

%%
alpha= VP;
alpha(S~=3)=0;
alpha= cube_alpha(alpha,size_pixel);
VP(S==3)=alpha(S==3);

%%
beta = VP;
beta(S~=5&S~=4)=0;
beta = cube_beta(beta,size_pixel);
VP(S==5|S==4)=beta (S==5|S==4);

%%
gamma = VP;
gamma(S~=6&S~=7)=0;
gamma = cube_gamma(gamma,size_pixel);
VP(S==6|S==7)=gamma(S==6|S==7); 
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
