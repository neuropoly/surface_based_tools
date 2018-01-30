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
% Function is written by C.Van Assel Univertité polytechnique de Montreal
% (November 2016)

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