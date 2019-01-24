# S3segmenter


# Dockerized app for segmentation of large tissue samples or TMA cores and output nuclei and/or cytoplasm masks with the option of exporting mean and median intensities of each segmented cell.

## Segmentation Methodology
This pipeline takes in a tiff stack of images. Specify the nuclei, cytoplasmic, and background channels in the file called S3segmenter_config.yml file. The UNet model has been trained on two classes : nuclei contours and everything else to increase the model's ability to identify nuclei contours for splitting. From the contour class probability map, a Laplacian of Gaussian (LoG) kernel is convolved over the image to identify regional maxima that approximate the center of each nuclei. These are used as seed points for a marker controlled watershed segmentation. For segmentating the cytoplasm, the nuclei are used as seed points for yet another marker controlled watershed segmentation. Although this produces some false positive objects in the background, these are eliminated based on size and intensity in the nuclei channel using robust measurements. The nuclei and cytoplasm label masks are saved as separate files along with mask edges overlayed on the original nuclei image.

Note: if segmenting TMA cores, it will be necessary to segment the core first to eliminate spurious objects and other cores that have crept into the image. This is done using coreSplitter.m that uses active contours and voronoi tesselation to identify the core of interest (usually near the center of the field of view). This pre-mask is used to clean up the nuclei and cytoplasm label masks later. In the config file, set preMask to 'true'. Set to 'false' if segmenting chunky tissue.

## Installation
Make sure Docker is installed (https://www.docker.com/products/docker-desktop)
Download the docker image using the command `docker pull clarenceyapp/s3segmenter' 

To build from source, download the code using `https://github.com/HMS-IDAC/S3segmenter.git`
Within Matlab, navigate to the S3segmenter/S3segmenterCODE folder.
The main wrapper function is in S3segmenterWrapper.m. To run this on one or over several .tif files, just use batchS3segmenterWrapper.m all the time. 
Using Matlab's compilation functionality, run the command `mcc -m S3segmenterWrapperDocker.m` to build the necessary binary
The binary can be run using Matlab Runtime (https://www.mathworks.com/products/compiler/matlab-runtime.html)
The random forest model ('model.mat') should exist in the dockerbuild/matlabDependencies folder. NOTE: THE RANDOM FOREST PORTION OF THE CODE DOES NOT WORK YET.
Alternatively, the code can be run by building a fresh docker container, using the Dockerfile contained in S3segmenter/dockerBuild/Dockerfile which installs all necessary dependencies
Then follow the usage instructions below for using the Docker container 

## Usage 
Create an input folder which contains any number of registered tiff stacks where channels are concatenated in the 3rd dimension.

Create a configuration folder, containing the file S3segmenter_config.yml
All configuration parameters must be set by the user and will be explained below later. Example configuration file can be found at 
https://github.com/HMS-IDAC/S3segmenter/dockerBuild/config/

Run the docker container, indicating paths to the input images, configuration file, and desired output location
`docker run -v /local/path/to/config/:/config -v /local/path/to/input/:/input -v /local/path/to/output/:/output clarenceyapp/s3segmenter`

## Input Data
Input data must be .tif image stacks. TMA cores or chunky tissues are acceptable. You must also include the preprocessed class probability map from the UNet model. This file should have the suffix '_NucSeg.tif'.

## Segmentation Output
Output folder must be specified, and can correspond to any already created local folder. Output files consist of subfolders - one for each .tif - containing the label masks for nuclei and cytoplasm. Also included are tab delimited files of median and mean intensities if Measure Features was selected in the config file.

## Configuration parameters
The following parameters must be set by the user based on their experimental design. 
Example configuration file can be found at: /localPath/S3segmenter/dockerBuild/config


Configuration parameters:

ip.addParamValue('ClassProbSource','unet',@(x)(ismember(x,{'RF','unet'})));
ip.addParamValue('NucMaskChan',[1],@(x)(numel(x) > 0 & all(x > 0 )));  the channel containing the nuclei marker
ip.addParamValue('CytoMaskChan',[2],@(x)(numel(x) > 0 & all(x > 0 )));  the channel containing a cytoplasm marker such as b-catenin
ip.addParamValue('TissueMaskChan',[3],@(x)(numel(x) > 0 & all(x > 0 )));  the channel containing a marker that generally stains the entire tissue. Could be autofluorescence even.
ip.addParamValue('preMask','true',@(x)(ismember(x,{'true','false'}))); % set to true if sample is TMA cores
ip.addParamValue('cytoMethod','distanceTransform',@(x)(ismember(x,{'RF','distanceTransform','ring'})));
ip.addParamValue('MedianIntensity','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('saveFig','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('saveMasks','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('segmentNucleus','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('measureFeatures','false',@(x)(ismember(x,{'true','false'}))); set to false to just save masks quickly and measure features in another program
ip.addParamValue('nucleiRegion','watershedContourInt',@(x)(ismember(x,{'watershedContourDist','watershedContourInt','dilation'})));
two methods of marker controlled watershed segmentation or just dilate around the nuclei.

ip.addParamValue('segmentCytoplasm','segmentCytoplasm',@(x)(ismember(x,{'segmentCytoplasm','loadMask','ignoreCytoplasm'})));
ip.addParamValue('useGPUArray','false',@(x)(ismember(x,{'true','false'}))); if you dont' have a gpu card, set to false.
ip.addParamValue('inferNucCenters','true',@(x)(ismember(x,{'true','false'}))); UNet produces a class probability map of contours and everything else. Setting inferNucCenters takes the inverse of the contour for later analysis.
ip.addParamValue('nucleiPriority','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('resizeFactor',1,@(x)(numel(x) == 1 & all(x > 0 )));  scale down to increase speed at the expense of losing image resolution
ip.addParamValue('logSigma',[2.5],@(x)(numel(x) >0 & all(x > 0 ))); the size of the LoG filter to approximate diameter of cell. Add multiple sizes if needed.
ip.addParamValue('upSample',2,@(x)(numel(x) == 1 & all(x > 0 )));  upsampling is needed to segment cytoplasm between closepacked nuclei.
ip.addParamValue('Docker','false',@(x)(ismember(x,{'true','false'}))); set to true if dockerizing
ip.addParamValue('dockerParams',0,@(x)(numel(x)==1)); leave this alone


