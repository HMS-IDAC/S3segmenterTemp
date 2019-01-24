function batchS3segmenterWrapper(testPath,varargin)


ip = inputParser;
ip.addParamValue('ClassProbSource','unet',@(x)(ismember(x,{'RF','unet'})));
ip.addParamValue('NucMaskChan',[1],@(x)(numel(x) > 0 & all(x > 0 )));  
ip.addParamValue('CytoMaskChan',[2],@(x)(numel(x) > 0 & all(x > 0 )));  
ip.addParamValue('TissueMaskChan',[3],@(x)(numel(x) > 0 & all(x > 0 )));  
ip.addParamValue('preMask','true',@(x)(ismember(x,{'true','false'}))); % set to true if sample is TMA cores
ip.addParamValue('cytoMethod','distanceTransform',@(x)(ismember(x,{'RF','distanceTransform','ring'})));
ip.addParamValue('MedianIntensity','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('saveFig','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('saveMasks','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('segmentNucleus','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('measureFeatures','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('nucleiRegion','watershedContourInt',@(x)(ismember(x,{'watershedContourDist','watershedContourInt','dilation'})));
ip.addParamValue('segmentCytoplasm','segmentCytoplasm',@(x)(ismember(x,{'segmentCytoplasm','loadMask','ignoreCytoplasm'})));
ip.addParamValue('useGPUArray','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('inferNucCenters','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('nucleiPriority','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('resizeFactor',1,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('logSigma',[2.5],@(x)(numel(x) >0 & all(x > 0 )));
ip.addParamValue('upSample',2,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('Docker','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('dockerParams',0,@(x)(numel(x)==1));
ip.parse(varargin{:});          
p = ip.Results;  

if isequal(p.Docker,'true')
    testPath = p.dockerParams.parentPath;
    p.modelPath = [p.dockerParams.modelPath filesep];
    p.modelCatPath = [p.dockerParams.modelCatPath filesep];
else
    if nargin<1 
        testPath = uigetdir(pwd,'Select a folder with image(s) to process');
    end
    p.modelPath ='';
    p.modelCatPath ='';
    p.outputPath = '';
end

%% read file names
    fileList = dir([testPath filesep '*.tif']);
    finalFileList = [];
    for iFile = 1:length(fileList)
        fName = fileList(iFile).name;
        if ~isfolder(fName) && ~contains(fName,'NucSeg')  && ~contains(fName,'._') ...
                && ~contains(fName,'.ome')
            finalFileList{end+1} = fName;
%                             
        end
    end
    disp (['Found ' num2str(length(finalFileList)) ' file(s)!'])
    
    for iFile = 1:length(finalFileList)
        S3segmenterWrapper(testPath,finalFileList{iFile},p);
    end
        