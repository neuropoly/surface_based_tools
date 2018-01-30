% developped by Gabriel Mangeat on 11/2016
% gabriel.mangeat@polymtl.ca
% Last modification: 16/07/2017

% the purpose of this code is to generate personalized atlaes of the
% cortical regions defined by a surfacic atlas in freesurfer space (here, 
% the PALS-B12 Brodmann area (BA) surface atlas (Van Essen 2005) was used)

clc 
clear all

fsdir = '/Applications/freesurfer/subjects/';
subjects_list = {'Ctrl028_R01_01' };
wanted_struct_PALS = {'???','Brodmann.20','Brodmann.21','Brodmann.22','Brodmann.27','Brodmann.28','Brodmann.35','Brodmann.36','Brodmann.37','Brodmann.38','MEDIAL.WALL','Brodmann.1','Brodmann.2','Brodmann.3','Brodmann.4','Brodmann.5','Brodmann.6','Brodmann.7','Brodmann.8','Brodmann.9','Brodmann.10','Brodmann.11','Brodmann.17','Brodmann.18','Brodmann.19','Brodmann.23','Brodmann.24','Brodmann.25','Brodmann.26','Brodmann.29','Brodmann.30','Brodmann.31','Brodmann.32','Brodmann.39','Brodmann.40','Brodmann.41','Brodmann.42','Brodmann.43','Brodmann.44','Brodmann.45','Brodmann.46','Brodmann.47'};
%wanted_struct_PALS = {'???','Brodmann.20','Brodmann.21','Brodmann.22','Brodmann.27','Brodmann.28','Brodmann.35','Brodmann.36','Brodmann.37','Brodmann.38','MEDIAL.WALL'};
hemi = {'rh','lh'};


% annotation to label:
%  mri_annotation2label  --annotation /Applications/freesurfer/subjects/fsaverage/label/lh.PALS_B12_Brodmann.annot
% --hemi lh --subject fsaverage --outdir /Applications/freesurfer/subjects/fsaverage/label/PALS_B12_BA_labels 
% idem rh ..


% mri_label2label
fsavg_labeldir = [fsdir 'fsaverage/label/PALS_B12_BA_labels'];

for iSub = 1:length(subjects_list)
    labeldir = [fsdir subjects_list{iSub} '/label/PALS_B12_BA_labels'];
    %unix(['mkdir ' labeldir]);
    
    for iBA = 1:length(wanted_struct_PALS)
        for iHemi = 1:length(hemi)
    
            cmd = ['mri_label2label --srcsubject fsaverage --srclabel ' ...
                fsavg_labeldir '/' hemi{iHemi} '.' wanted_struct_PALS{iBA} '.label ' ...
                 '--trgsubject ' subjects_list{iSub} ' --hemi ' hemi{iHemi} ...
                 ' --trglabel ' labeldir '/' hemi{iHemi} '.' wanted_struct_PALS{iBA} '.label ' ...
                 '--regmethod surface'];
            disp(cmd);
            %unix(cmd);
            
        end
    end
end

% %%%%%%%%%%%%
% % mris_label2annot
% ctabfile = '/Applications/freesurfer/subjects/fsaverage/label/PALS_B12.annot.ctab';
% 
% 
% % for iHemi = 1:length(hemi)
% % 
% %     for iSub = 1:length(subjects_list)
% %         labeldir = [fsdir subjects_list{iSub} '/label/PALS_B12_BA_labels'];
% %         
% %         label_list = [];
% %         for iBA = 1:length(wanted_struct_PALS)
% %             label_list = [label_list '--l  ' labeldir '/' hemi{iHemi} '.' wanted_struct_PALS{iBA} '.label '];
% %         end
% %         
% %         cmd = ['mris_label2annot --s ' subjects_list{iSub} ...
% %            ' --h ' hemi{iHemi} ...
% %            ' --ctab ' ctabfile ...
% %            ' ' label_list ...
% %            ' --a ' hemi{iHemi} 'PALS_B12_Brodmann --maxstatwinner'];
% %         unix(cmd);
% %         %disp(cmd);
% %         
% %     end 
% % end
% %%%%%%%%%%%%%%% 

% mri_label2vol 
path_r01 = '/Users/gabriel/Desktop/r01/';
template_name = 'ep2d_diff_qb64_b1k_eddy.nii.gz';
regmat_name = 'diff_first_b0_disco_2_fsT1_12dofbbr.dat';


for iSub = 1:length(subjects_list)
    labeldir = [fsdir subjects_list{iSub} '/label/PALS_B12_BA_labels'];

    unix(['mkdir ' path_r01 subjects_list{iSub} '/atlas']);
    unix(['mkdir ' path_r01 subjects_list{iSub} '/atlas/sreg_BAs']);

    for iBA = 1:length(wanted_struct_PALS)

        for iHemi = 1:length(hemi)
            cmd = ['mri_label2vol --label ' labeldir '/' hemi{iHemi} '.' wanted_struct_PALS{iBA} '.label '...
                '--temp ' path_r01 subjects_list{iSub} '/diff/' template_name ...
                ' --reg ' path_r01 subjects_list{iSub} '/mat/' regmat_name ...
                ' --fillthresh 1 --proj frac -0.1 1 0.1 ' ...
                ' --o ' path_r01 subjects_list{iSub} '/atlas/sreg_BAs/' hemi{iHemi} '_'  wanted_struct_PALS{iBA} '.nii.gz'  ...
                ' --subject ' subjects_list{iSub} ' --hemi ' hemi{iHemi}];
            disp(cmd);
            %unix(cmd);
            
        end
    end
end 

%%% END of generating BA atlases %%%
%%% generating aseg atlases %%%


% sending aseg in diff space
aseg_file = 'aseg.auto_noCCseg';

for iSub = 1:length(subjects_list)
    cmd = ['mri_label2vol --seg ' fsdir subjects_list{iSub} '/mri/' aseg_file '.mgz' ...
        ' --temp ' path_r01 subjects_list{iSub} '/diff/' template_name ...
        ' --reg ' path_r01 subjects_list{iSub} '/mat/' regmat_name ...
        ' --o ' fsdir subjects_list{iSub} '/mri/aseg_file_DIFF_space.nii.gz'];
    disp(cmd);
    unix(cmd);
        
end


% extracting roi of interest from aseg
wanted_struct_aseg = {'Left_Cerebellum_White_Matter' 'Left_Thalamus' 'Left_Caudate' 'Left_Putamen' 'Brain_Stem' 'Left_Hippocampus' 'Left_Amygdala'...
                     'Right_Cerebellum_White_Matter' 'Right_Thalamus' 'Right_Caudate' 'Right_Putamen' 'Right_Hippocampus' 'Right_Amygdala'};
index_struct_aseg = [7 9 11 12 16 17 18 46 48 50 51 53 54];
prefix = 'struct_';


for iSub = 1:length(subjects_list)
    for iStruct = 1:length(wanted_struct_aseg)
        cmd = ['fslmaths ' fsdir subjects_list{iSub} '/mri/aseg_file_DIFF_space.nii.gz' ...
            ' -uthr ' num2str(index_struct_aseg(iStruct)) ' -thr ' num2str(index_struct_aseg(iStruct)) ...
            ' ' path_r01 subjects_list{iSub} '/atlas/sreg_BAs/' prefix wanted_struct_aseg{iStruct} '.nii.gz'];
        disp(cmd);
        unix(cmd);  
        
        cmd = ['fslmaths ' path_r01 subjects_list{iSub} '/atlas/sreg_BAs/' prefix wanted_struct_aseg{iStruct} '.nii.gz' ...
            ' -div ' num2str(index_struct_aseg(iStruct)) ...
            ' ' path_r01 subjects_list{iSub} '/atlas/sreg_BAs/' prefix wanted_struct_aseg{iStruct} '.nii.gz'];
        disp(cmd);
        unix(cmd);  
        
    end
end



% fslmaths /[path to aseg file]/s001_aseg2raw.nii.gz -uthr 11 -thr 11 /[path to working directory]/ROIS/l_caudate_freesurfer_rawavg.nii.gz

% mri_label2vol --seg /Applications/freesurfer/subjects/Ctrl028_R01_01/mri/aseg.auto_noCCseg.mgz --temp /Users/gabriel/Desktop/r01/Ctrl028_R01_01/diff/ep2d_diff_qb64_b1k_eddy.nii.gz --reg /Users/gabriel/Desktop/r01/Ctrl028_R01_01/mat/diff_first_b0_disco_2_fsT1_12dofbbr.dat --o /Users/gabriel/Desktop/r01/Ctrl028_R01_01/test.nii

% mri_concat % if need to put back every roi in one nii file

