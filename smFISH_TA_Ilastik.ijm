// Created by Yew Yan Wong 23 November 2022
dir1 = getDirectory("Choose Source Directory ");
dir2 = getDirectory("Choose Destination Directory");
list = getFileList(dir1);
bg_value = 500; // This is the vValue for background correction. This value (500) needs to be changed depending on the microscope used. ItThis value should be aboutalmost the same as the average background fluorescence intensity of a place where no cells and debris are contained, after the imaging conditions are fixed.
pixel_x = 0.13; // Enter the length in the x-axis per pixel.
pixel_y = 0.13; // Enter the length in the y-axis per pixel.
for (i=0; i<list.length; i++) {
  showProgress(i+1, list.length);
  run("Bio-Formats", "open="+dir1+"/"+list[i]+" color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");rename("stack"); // Importing Images.
  run("Properties...", "channels=3 slices=1 frames=10 pixel_width="+pixel_x+" pixel_height="+pixel_x+" voxel_depth=1.0000000 frame=[0.00 sec]"); // Change the image properties. The pixel width setting and the order of the images are changed. The acquired image stack is a three-color time-lapse image, but here it is recognized as a three-color z-stack image. This enables the average projection of the image at athe later step.
  run("Gaussian Blur...", "sigma=1 stack"); // Apply a Ggaussian blur filter.
  run("Split Channels"); // Separate the image stack into channels. 
  selectWindow("C1-stack"); // Select an image for the first color (SNAPtag channel).
  run("Bleach Correction", "correction=[Simple Ratio] background="+ bg_value +""); // Correct photobleaching.
  run("Subtract Background...", "rolling=50 stack"); //  Correction of uneven illumination.
  selectWindow("C1-stack"); close(); // Delete the pre-processed image.
  selectWindow("C2-stack"); // Select an image for the second channel (RFP channel).
  run("Bleach Correction", "correction=[Simple Ratio] background="+ bg_value +""); // Correct photobleaching.
  run("Subtract Background...", "rolling=50 stack"); // Correct for uneven illumination.
  run("Subtract Background...", "rolling=5 stack"); // Emphasize spots.
  selectWindow("C2-stack"); close(); // Delete the pre-processed image.
  selectWindow("C3-stack"); // Select an image for the third channel (mNG channel). 
  run("Bleach Correction", "correction=[Simple Ratio] background="+ bg_value +""); // Correct photobleaching.
  run("Subtract Background...", "rolling=50 stack"); // Correct for uneven illumination. 
  run("Subtract Background...", "rolling=5 stack"); // Emphasize spots.
  selectWindow("C3-stack"); close(); // Delete the pre-processed image.
  run("Merge Channels...", "c1=DUP_C1-stack c2=DUP_C2-stack c3=DUP_C3-stack create ignore"); // Merge the processed three- color stacks.
  run("Z Project...", "projection=[Average Intensity]"); // Average intensity projection
  close("Composite"); close("Merged"); // Delete the unwanted images.
  // rename("Nst-SNAP-Brd4"); // Rename the image.
  run("Split Channels");
  saveAs("TIFF", dir2+substring(list[i], 0, list[i].length - 4) + "_mNG");
  close();
  saveAs("TIFF", dir2+substring(list[i], 0, list[i].length - 4) + "_RFP");
  close();
  saveAs("TIFF", dir2+substring(list[i], 0, list[i].length - 4) + "_SNAP");
  close();
}