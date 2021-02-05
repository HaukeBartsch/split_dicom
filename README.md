# Split DICOM images into left and right

Animal imaging routinely captures images from more than one individual. This is done mostly to shorten the acquisition time, utilize PET trackers better. In order to map one animal to one DICOM study this project receives a DICOM image series and generates two new image series in the output directory of the left and rigth half of each 2D DICOM slice (assuming an axial oriented slice aquisition.