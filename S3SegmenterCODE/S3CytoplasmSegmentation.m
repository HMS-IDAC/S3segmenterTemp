function [cytoplasmMask,nucleiMask] = S3CytoplasmSegmentation(nucleiMask,cyto,tissue,modelCat,varargin)

ip = inputParser;
ip.addParamValue('cytoMethod','distanceTransform',@(x)(ismember(x,{'RF','distanceTransform','ring'})));
ip.addParamValue('useGPUArray','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('nucleiPriority','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('resize',1,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('sizeFilter',1,@(x)(numel(x) == 1 & all(x > 0 ))); 
ip.addParamValue('upSample',2,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('mask', [], @(x) isnumeric(x) || islogical(x));
ip.parse(varargin{:});          
p = ip.Results;  
        
%% mask
if ~isempty(p.mask)
    cyto = cyto.*imresize(cast(p.mask,class(cyto)),size(cyto));
    tissue = tissue.*imresize(cast(p.mask,class(tissue)),size(tissue));       
end

%% cytoplasm segmentation methods
        
        switch p.cytoMethod
            case 'RF'
          F = pcImageFeatures(imresize(double(cyto)/65335,p.resize,'bilinear'),modelCat.sigmas,modelCat.offsets,...
              modelCat.osSigma,modelCat.radii,modelCat.cfSigma,modelCat.logSigmas,modelCat.sfSigmas,...
              modelCat.ridgeSigmas,modelCat.ridgenangs,modelCat.edgeSigmas,modelCat.edgenangs,modelCat.nhoodEntropy,...
              modelCat.nhoodStd);
             [imL,catClassProbs] = imClassify(F,modelCat.treeBag,100);
              contours = imresize(catClassProbs(:,:,2),2);
              bgm =imresize(bwmorph( imgaussfilt3(cyto,2)<100,'thin',Inf),2);
              cytograd= imimposemin(imresize(contours,[size(nucleiMask,1) size(nucleiMask,2)]),bgm|nucleiMask);
              cellMask= watershed(cytograd);

            case 'distanceTransform'
                tissueProc=imgaussfilt3(imresize(tissue,0.25),1);
                NeighborhoodSize  =round(2*sqrt(p.sizeFilter/pi))*25;
                if mod(NeighborhoodSize  ,2) == 0
                    NeighborhoodSize=NeighborhoodSize+1;
                end
                Th= adaptthresh(tissueProc,0.8,'NeighborhoodSize',NeighborhoodSize);
                preCellMask = (imresize(imbinarize(tissueProc,Th),size(nucleiMask),'nearest')+nucleiMask)>0;
                gdist = graydist(cyto,nucleiMask);
                cytograd = imimposemin(gdist, (imerode(1-(preCellMask+nucleiMask)>0,strel('disk',5))>0) | nucleiMask );
                cellMask=watershed(cytograd);
%                 cellMask=cellMask.*cast(preCellMask,class(cellMask));
%                 cellMask = bwlabel((nucleiMask + bwareaopen(cellMask,50))>0);
              

%                  nMaskDist =-bwdist(~nucleiMask);
%                  cytograd= imimposemin(imresize(nMaskDist,[size(nucleiMask,1) size(nucleiMask,2)]),nucleiMask);
%                  cellMask=watershed(cytograd);

            case 'contours'
                
                contours = normalize(steerableDetector(im2double(cyto),2,1.5));

            case 'ring'
                cellMask = bwlabel(imdilate(nucleiMask,strel('square',3)));
        end


      

        %% eliminate 'cells' without nuclei
        if isequal(p.nucleiPriority,'true' )
            %clean this up!
            test=cast(~ismember(unique(cellMask),unique(nucleiMask)),class(unique(cellMask))).*unique(cellMask);
            bgCells=find(test>0);
            for iTest=bgCells'
                cellMask(cellMask ==test(iTest))=0;
            end
        end


        %% eliminate border cells
%         inCells = imclearborder(cellMask>0);
%         borderCells = cellMask.*cast((cellMask>0)-inCells,class(cellMask));
%         borderIdx = unique(borderCells);
%         nucleiMask_border = ~ismember(bwlabel(nucleiMask),borderIdx);
% 
%         cellMask = bwlabel(inCells);
%         nucleiMask = nucleiMask_border.*cellMask;
%         cytoplasmMask = cellMask - nucleiMask;
%         
        %% filter based on a nuclei
          stats=regionprops(cellMask,nucleiMask>0,'MaxIntensity','Area');
          idx = find([stats.MaxIntensity] > 0 );
          finalCellMask = bwlabel(ismember(cellMask,idx));
          nucleiMask = cast(nucleiMask,class(finalCellMask)).*finalCellMask; 
          cytoplasmMask = finalCellMask - nucleiMask;
        