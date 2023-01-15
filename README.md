# kidney-segmentation

Description:
This Matlab code provides a graphical user interface (GUI) for segmenting the kidneys in a semi-automatic manner. The GUI allows the user to load a set of kidney MR image, selects both the kidneys semi-automatically, and saves the segmented image as binary masks. 

Requirements:
Matlab 2016b or later; 
Image processing toolbox; 
MRELab (for further MRE data processing)

Usage:
Run the "MRE_KIDNEY_DEEP.m" file in Matlab to open the GUI. Use the "Load Images" button to select a set of kidney slices (magnitude and phase kidney MRE DICOM images). Seletc the magnitude kidney MRE images. Use the "Select Region" button to select the kidneys in the image. You can use the "add" and "remove" options to change the selected regions. You can also select the "biggest" or the "smallest" segmented regions. You can also use the "Reset" button to clear the current selection and start over.

Input:
A series of kidney MRE magnitude and phase images

Output:
Kidney binary masks and vars files

Authors:
Deep B. Gandhi
