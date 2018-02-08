function cmd = pvc_create_expanded_surfaces( surf, fsdir, N, voxel_size)

% This function create 2*N expanded surfaces, N in each side of the
% original surface. 
%
% Inputs: 
%           surf: name of the input surface
%           N: number of expanded surfaces to be created
%           fsdir: directory that contains the surface
%           voxel_size: size 
%
% Outputs:
%           2*N expanded surfaces are created in the fsdir directory
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

current_dir = pwd;

cd(fsdir);
launch_bash_profile = 'source ~/.bash_profile;';

cmd = launch_bash_profile;
for iter=1:N 
    % We create N intermediary surfaces
    % The max distance from the first surface to the last expanded one has 
    % to be sqrt(3), which is the maximale diagonale of a pixel
    exp = iter*voxel_size*1.6/N; 
    if iter<10
        cmd = horzcat(cmd, ['mris_expand ' surf ' ' num2str(exp) ' ' surf '00' num2str(iter) ';']);  
    else
        cmd = horzcat(cmd, ['mris_expand ' surf ' ' num2str(exp) ' ' surf '0' num2str(iter) ';']); 
    end
end 
for iter=1:N 
    % We create n intermediary surfaces in the other side
    exp = -iter*voxel_size*1.6/N;
    if iter<10
        cmd = horzcat(cmd, ['mris_expand ' surf ' ' num2str(exp) ' ' surf 'n00' num2str(iter) ';']);  
    else
        cmd = horzcat(cmd, ['mris_expand ' surf ' ' num2str(exp) ' ' surf 'n0' num2str(iter) ';']); 
    end
end
disp(cmd);
unix(cmd);

cmd = ['cp ' surf ' ' surf '000'];
disp(cmd);
unix(cmd);

cd(current_dir);

end