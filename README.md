# Split DICOM images into left and right

Animal imaging routinely captures images from more than one individual. This is done mostly to shorten the acquisition time and utilize PET trackers better. In order to map one animal to one DICOM study this project receives a DICOM image series and generates two new image series in the output directory for the left and rigth half of each 2D DICOM slice (assuming an axial oriented slice aquisition).

## Build

There are docker build instructions available:

```
git clone https://github.com/HaukeBartsch/split_dicom.git
cd split_dicom
docker build -t split_dicom build
```
## Run

For development its nicer to mount the source folder inside the docker container. The build can be done in Linux which the hosts editor can be used to develop features:

```
docker run --rm -it -v /Users/hauke/src/split_dicom:/root/split_dicom split_dicom /bin/bash
cd /root/split_dicom
cmake .
make
./split_dicom
```