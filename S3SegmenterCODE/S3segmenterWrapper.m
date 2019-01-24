function S3segmenterWrapper(testPath,fileName,p)
%this function requires input of nuclei stack range. It assumes that every
%stack beyond that to the end is a cytoplasmic stain. Marker controlled
%watershed based on distance transform of nuclei channel is employed to
%separate nuclei clumps.

% training with nuclei model = pixelClassifierTrain('Y:\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\20x full exp\nucleiTrainingSet\halfsize','sigmas',[3 7 15 30],'sfSigmas',[],'edgeSigmas',1,'ridgeSigmas',1,'radii',[3 5 10 ],'logSigmas',[4 6 8],'nhoodStd',[5 11 15],'pctMaxNPixelsPerLabel',100,'adjustContrast',0);
% training with nuclei model = pixelClassifierTrain('Y:\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\Real experiment\All Cyles\nucleiTrainingSet10x','sigmas',[1 3 7 9],'sfSigmas',[],'edgeSigmas',2,'ridgeSigmas',2,'radii',[3 5 7 ],'logSigmas',[2 4 6 ],'nhoodStd',[5 11 15],'pctMaxNPixelsPerLabel',100,'adjustContrast',0);

% training with nuclei modelCat = pixelClassifierTrain('Y:\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\20x full exp\cateninTrainingSet\halfsize','sigmas',[1 3 7 ],'sfSigmas',[1 2 3],'logSigmas',[2 4],'nhoodStd',[3 5 11],'pctMaxNPixelsPerLabel',100,'adjustContrast',0);
% training with nuclei modelCat = pixelClassifierTrain('Y:\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\20x full exp\cateninTrainingSet\RFTrainingSet','sigmas',[1 3 5 7  ],'sfSigmas',[],'edgeSigmas',[2],'ridgeSigmas',[1 2],'logSigmas',[1 2 4 ],'nhoodStd',[5 7 11 ],'pctMaxNPixelsPerLabel',40,'adjustContrast',0);
% training with nuclei modelCatSF1 = pixelClassifierTrain('Y:\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\20x full exp\cateninTrainingSet','sigmas',[1 3 7 ],'sfSigmas',[1],'logSigmas',[1 2 4 6],'nhoodStd',[5 11 15],'pctMaxNPixelsPerLabel',100,'adjustContrast',0);
% training with nuclei modelCat = pixelClassifierTrain('Y:\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\20x full exp\cateninTrainingSet\RFTrainingSet','sigmas',[1 3 5 7 ],'sfSigmas',[1 2],'logSigmas',[1 2 4 6],'nhoodStd',[5 11 15],'pctMaxNPixelsPerLabel',50,'adjustContrast',0);
% training with nuclei modelCat = pixelClassifierTrain('Y:\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\20x full exp\cateninTrainingSet\manualRFTrainingSet','sigmas',[1 2 4 ],'sfSigmas',[1 2],'edgeSigmas',[1 2],'ridgeSigmas',[1 2],'logSigmas',[1 2 4],'nhoodStd',[3 5 11 15],'pctMaxNPixelsPerLabel',100,'adjustContrast',0);


% override default values because using Docker config file
modelPathName = p.modelPath;

if nargin < 1 
    if nargin <2 
         testPath = pwd;
    end
     [fileName, testPath] = uigetfile([testPath filesep '*.tif'],'Select file to process');
else
     testPath = [testPath filesep];
end
[~,filePrefix] = fileparts(fileName);


    %% generate class probability maps
    switch p.ClassProbSource
        case 'RF'
            if isempty(modelPathName)
                [modelFileName, modelPathName] = uigetfile([testPath filesep '*.mat'],'Select RF model');
            end
% modelFileName = 'modelReducedFeatures.mat';
% modelPathName = 'D:\LSP\cycif\S3\nucleiTrainingSet\'; 
        %#pragma treeBagger
        load([modelPathName modelFileName])

        nuclei = imread([testPath fileName]);
        [nucleiPM,nucleiCrop]=RFSplitter(nuclei,modelNuc,'split',false);
        case 'unet'
        
        nucleiCrop = imread([testPath fileName],p.NucMaskChan);
        nucleiPM = imresize(volumeRead([testPath filePrefix '_NucSeg.tif']),p.resizeFactor);
        nucleiPM = imresize(nucleiPM(:,:,2),size(nucleiCrop));
    end
      
tic

%% mask the core
if isequal(p.preMask,'true')
    TMAmask = coreSegmenter(nucleiCrop);
else
    TMAmask = [];
end

   %% nuclei segmentation
   if isequal(p.segmentNucleus,'true')
%        nucleiPM=S3tileReturn(nucleiPM);
%        nucleiCrop = S3tileReturn(nucleiCrop);
       
        [nucleiMask,largestNucleiArea] = S3NucleiSegmentation(nucleiPM,nucleiCrop,p.logSigma,'mask',TMAmask,...
            'inferNucCenters',p.inferNucCenters,'nucleiRegion',p.nucleiRegion,'resize',p.resizeFactor,'upSample',p.upSample);
   else
        nucleiMask = imread([testPath filesep name '_nucleiLM' ext]);
   end


    %% cytoplasm segmentation
   switch p.segmentCytoplasm
       case 'segmentCytoplasm'
            cyto =[];
            for iChan = p.CytoMaskChan
               cyto= cat(3,cyto,imread([testPath fileName],iChan));
            end
            cyto = max(cyto,[],3);
            % load random forest model if 
            if isequal(p.cytoMethod,'RF')
                %#pragma treeBagger
                load('D:\LSP\cycif\S3\cytoTrainingSet\modelContours1.mat')
            else
            modelCat=[];
            end
            
%             cyto=S3tileReturn(cyto);
            cyto = imresize(cyto,size(nucleiMask));
            tissue = imread([testPath fileName],p.TissueMaskChan);
%             tissue = S3tileReturn(tissue);
            tissue = imresize(tissue,size(nucleiMask));
            [cytoplasmMask,nucleiMask]=S3CytoplasmSegmentation(nucleiMask,cyto,tissue,modelCat,'mask',TMAmask,...
                'cytoMethod',p.cytoMethod,'resize',p.resizeFactor,'sizeFilter',largestNucleiArea);

        case 'loadMask'
            nucleiMask = bwlabel(nucleiMask);
            cytoplasmMask = imread([testPath filesep '_' name '_cytoLM' ext]);
        case 'ignoreCytoplasm'
            nucleiMask = bwlabel(nucleiMask);
            cytoplasmMask = nucleiMask;
   end
   
   %% set up output directories
       
if isequal(p.Docker,'true')
    outputPath = p.dockerParams.outputPath;
else
    if exist([testPath  'output'],'dir')~=7
        mkdir([testPath  'output'])
    end
    outputPath = [testPath 'output'];
end

[~,name]=fileparts(fileName);
if exist([outputPath filesep name],'dir')~=7
    mkdir([outputPath filesep name])
end
outputPath = [outputPath filesep name];  

    %% measureFeatures
    if isequal(p.measureFeatures,'true')
        S3MeasureFeatures(nucleiMask,cytoplasmMask,testPath,fileName,'MedianIntensity',p.MedianIntensity,'outputPath',outputPath,'Docker',p.Docker);
    end

   
    %% display
   
    if isequal(p.saveFig,'true')
        if ~isequal(p.Docker,'true')
            imshow(imresize(nucleiCrop,[size(nucleiMask,1) size(nucleiMask,2)]),[]), hold on, visboundaries(bwboundaries(nucleiMask),'LineWidth',1)
            savefig ([outputPath filesep name '_nucleiMasked.fig' ])
        end
%         savefig ([analysisPath filesep name '_nucleiMasked.fig' ])
        % add mask index to image

        % save image outlines for overlaying in OMERO
        tiffwriteimj(cat(3,uint16(bwperim(nucleiMask))*max(nucleiCrop(:)),imresize(nucleiCrop,[size(nucleiMask,1) size(nucleiMask,2)],'nearest')),[outputPath filesep name '_nucleiOutlines.tif'])
        tiffwriteimj(cat(3,uint16(bwperim(cytoplasmMask))*max(nucleiCrop(:)),imresize(nucleiCrop,[size(cytoplasmMask,1) size(cytoplasmMask,2)],'nearest')),[outputPath filesep name '_cytoplasmOutlines.tif'])                            
        tiffwriteimj(nucleiMask,[outputPath filesep name '_nucleiMask.tif'])
        tiffwriteimj(cytoplasmMask,[outputPath filesep name '_cytoplasmMask.tif'])

        close all
        disp(['Completed ' fileName])
    end
    toc 
end


