function [nucleiMask,largestNucleusArea] = S3NucleiSegmentation(NucleiPM,nucleiImage,logSigma,varargin)

ip = inputParser;
ip.addParamValue('nucleiRegion','watershedContourInt',@(x)(ismember(x,{'watershedContourDist','watershedContourInt','dilation'})));
ip.addParamValue('useGPUArray','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('inferNucCenters','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('resize',1,@(x)(numel(x) == 1 & all(x > 0 ))); 
ip.addParamValue('upSample',2,@(x)(numel(x) == 1 & all(x > 0 ))); 
ip.addParamValue('mask', [], @(x) isnumeric(x) || islogical(x));
ip.parse(varargin{:});          
p = ip.Results;  

%% mask
if ~isempty(p.mask)
    nucleiImage = nucleiImage.*imresize(cast(p.mask,class(nucleiImage)),size(nucleiImage));
    NucleiPM = NucleiPM.*imresize(cast(p.mask,class(NucleiPM)),size(NucleiPM));
end


%% preprocess
if size(NucleiPM,3)==2
    nucleiCentersResized = imresize(NucleiPM(:,:,2),p.resize);
else
    nucleiCentersResized = [];
end
nucleiContoursResized = imresize(NucleiPM(:,:,1),p.resize);
nucleiImageResized = imresize(nucleiImage,p.resize);



%% use nuclei contour class prob maps from UNet or RF
    if isequal(p.inferNucCenters,'false') || ~isempty(nucleiCentersResized) %RF
        nucleiClass = nucleiCentersResized;
    else % UNet
        nucleiClass = 1-nucleiContoursResized;
    end

    %% markers based on log filter on classProbs 3
     if isequal(p.useGPUArray,'true')
      logNuclei=  gather(imgaussfilt3(gpuArray(filterLoG(nucleiClass,logSigma)),2)); %3.5, 2
     else
         
      fgFiltered=[];
        for iLog = 1:numel(logSigma)
            fgFiltered(:,:,iLog) = filterLoG(nucleiClass,logSigma(iLog));
        end
     logNuclei= imhmax(sum(fgFiltered,3),1.5); %3.5, 2
     end
     logfgm = imregionalmax(logNuclei);
%      imshowpair(logfgm,nucleiImage)

    %% upsample if needed
    if p.upSample > 1
     nucleiContoursResized = imresize(nucleiContoursResized,p.upSample);
    end

%% apply watershed transform
switch p.nucleiRegion
    case 'watershedContourInt'
      gradmag2= imimposemin(nucleiContoursResized,imresize(logfgm,[size(nucleiContoursResized,1) size(nucleiContoursResized,2)],'nearest'));
      foregroundMask= watershed(gradmag2);

    case 'dilation'
      foregroundMask = imdilate(logfgm,strel('square',2));

    case 'watershedContourDist'
        
        gdist = graydist(nucleiContoursResized,imresize(logfgm,2));
        cytograd= imimposemin(gdist,imresize(logfgm,size(gdist)));
        foregroundMask=watershed(cytograd);
        
        
%         logDist = -bwdist(imresize(logfgm,[size(nucleiContoursResized,1) size(nucleiContoursResized,2)],'nearest'));
%         cytograd= imimposemin(logDist,imresize(logfgm,size(logDist)));
% 
%         foregroundMask=watershed(cytograd);
%         nucleiBlur=imgaussfilt3(nucleiImage,2);
%         foreground =nucleiBlur>thresholdOtsu(nucleiBlur);
%         foregroundMask = foregroundMask.*cast(imresize(foreground,[size(foregroundMask,1) size(foregroundMask,2)]),class(foregroundMask));

 end

    %% process mask
   allNuclei = foregroundMask;

   if isequal(p.inferNucCenters,'false')
       nucleiCentersResized = imresize(nucleiCentersResized,size(nucleiContoursResized));
    stats=regionprops(bwlabel(allNuclei),nucleiCentersResized+nucleiContoursResized,'MeanIntensity','MinIntensity','Area');   
    idx = find([stats.MeanIntensity] > 0.75 & [stats.Area] > 10 & [stats.Area] < median(cat(1,stats.Area))*5 );
   else
    stats=regionprops(bwlabel(allNuclei),imresize(nucleiImageResized,[size(allNuclei,1) size(allNuclei,2)]),'MeanIntensity','Area');
    MITh = median(cat(1,stats.MeanIntensity))-1.4286*0.5*mad(cat(1,stats.MeanIntensity));
    idx = find([stats.MeanIntensity] > MITh & [stats.Area] > 10 ...
        & [stats.Area] < median(cat(1,stats.Area))*5 );
    end

   nucleiMask = ismember(bwlabel(allNuclei),idx);
   statsNM=regionprops(imresize(nucleiMask,[size(allNuclei,1) size(allNuclei,2)]),'Area');
   largestNucleusArea=prctile([statsNM.Area],95);


