# partial_volume_correction
Partial Volume Effect (PVE) hampers the accuracy of studies aiming at mapping MRI signal in the cortex due to the close proximity of adjacent white matter (WM) and cerebrospinal fluid (CSF). The proposed framework addresses this issue by disentangling the various sources of MRI signal within each voxel, assuming three classes (WM, gray matter, CSF) within a small neighbourhood. This tool allow accurate extraction of MRI metrics using surface-based analysis. This tool can be particularly useful for probing pathology in outer or inner cortical layers, which are subject to strong PVE with adjacent CSF or WM.

# Dependencies
* freesurfer 5.3
* matlab R2016_a

# Installation
* Download the above functions
* Add the folder code to your Matlab path
* Launch a matlab tab from the console
* Run the function *pvc_partial_volume_estimation*

# Inputs
* MRI volume 
* Subject directory path

__OPTIONS__:
* reg file (if the volume has been registered on the surface)
* path for output

__METHOD OPTION__:
* second correction: default false. Can only be added if the volume is isotropic. Correspond to the non linearity of the volume fraction depending on the orientation of the surface normal(for additionnal informations, see the abstract in the DOC directory).

# Outputs
* White matter mask
* Grey matter mask
* CSF mask
* r2 map in specified path or in current folder

# Examples
If the mri volume is the same resolution than the aseg file:
pvc_partial_volume_estimation(‘/sujet_test/mri/T1.mgz’, ‘/sujet_test’, [1 1 1]);

If the resolution of the mri volume is different from the aseg file:
pvc_partial_volume_estimation(‘/sujet_test/mri/T1.mgz’, ‘/sujet_test’, [1 1 1], ‘subject_test/mri/aseg.mgz, ‘subject_test/mri/T1_surfreg.reg’);


