function TMAmask = coreSegmenter(DAPI)
% using active contours after entropy filtering
% for iFile = 21:122
%     DAPI = imread([num2str(iFile) '_NucSeg.tif'],1);
        nucGF=entropyfilt(imresize(DAPI,[100 100]),ones(5,5));
        buffer = round(0.2*size(nucGF,1));
        mask = zeros(size(nucGF));
        mask(buffer:end-buffer,buffer:end-buffer) = 1;
        nuclearMask = activecontour(nucGF,mask,'Chan-Vese','SmoothFactor',5);

        %% watershed segmentation
        nMaskDist =imgaussfilt3(-bwdist(~nuclearMask),2);
        cytograd= imimposemin(nMaskDist,imerode(~nuclearMask,strel('disk',3))| imregionalmin(nMaskDist));
        TMAmask=imclearborder(watershed(cytograd)>0);
        TMAlabel = bwlabel(TMAmask);
        stats= regionprops(TMAlabel);
        
        %% choose the object closest to center of image
        for iObject = 1: numel(stats)
            dist(iObject) = sqrt((stats(iObject).Centroid(1)-50)^2+(stats(iObject).Centroid(2)-50)^2);
        end
        [minDistance, indexMin] = min(dist);
        TMAmask = imdilate(imfill(imclose(TMAlabel==indexMin,strel('disk',5)),'holes'),strel('disk',3));
%        close all
%         imshowpair(imresize(DAPI,size(nucGF)),bwperim(TMAmask))
%         pause(0.1)
end