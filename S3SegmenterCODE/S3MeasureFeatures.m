function S3MeasureFeatures(nucleiMask,cytoplasmMask,testPath,fileName,varargin)   
ip = inputParser;
ip.addParamValue('MedianIntensity','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('Docker','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('outputPath',[],@isstr);
ip.parse(varargin{:});          
p = ip.Results; 
   
cytoChanEnd = numel(finalChannelList);
  
if isempty(p.outputPath)
    [~,name]=fileparts(fileName);
    analysisPath = [testPath 'output' filesep name];
    mkdir(analysisPath)
end
   

%% measure intensities from regions
            if max(nucleiMask(:))>max(cytoplasmMask(:))
                numCells = max(nucleiMask(:));
            else
                numCells = max(cytoplasmMask(:));
            end
            meanIntNucTable = zeros(numCells,cytoChanEnd);
            meanIntCytoTable = zeros(numCells,cytoChanEnd);
            medianIntNucTable = zeros(numCells,cytoChanEnd);
            medianIntCytoTable = zeros(numCells,cytoChanEnd);


            for iChan = 1: cytoChanEnd
                I= imresize(imread([testPath fileName],iChan ),[size(nucleiMask,1) size(nucleiMask,2)]);
                nucleiStats=regionprops(nucleiMask,I,'MeanIntensity','Centroid','Area','PixelIdxList');
                cytoStats=regionprops(cytoplasmMask,I,'MeanIntensity','Centroid','Area','PixelIdxList');

                if isequal(p.MedianIntensity,'true')

                    for iCell = 1: numel(nucleiStats)
                        medianIntNucTable(iCell,iChan) = median(I(nucleiStats(iCell).PixelIdxList));
                        medianIntCytoTable(iCell,iChan) = median(I(cytoStats(iCell).PixelIdxList));
                    end
                end

                meanIntNucTable(:,iChan) = [nucleiStats.MeanIntensity]';
                meanIntCytoTable(:,iChan) = [cytoStats.MeanIntensity]';
                disp(['Measured channel ' int2str(iChan)])
            end

            meanIntTable = [meanIntNucTable meanIntCytoTable];
            medianIntTable = [medianIntNucTable medianIntCytoTable];
            areaTable = [cat(1,nucleiStats.Area) cat(1,cytoStats.Area)  ];
            centroidCellTable = cat(1,nucleiStats.Centroid);       

            %% write results to txt file
            if ~isempty(areaTable)
                variableNucNamesMeanIntensity = {};
                variableCytoNamesMeanIntensity = {};
                variableNucNamesMedianIntensity = {};
                variableCytoNamesMedianIntensity = {};

                if exist([testPath 'channel_metadata.csv']) ==2
                    [~,~,channelNames] = xlsread([testPath 'channel_metadata.csv']);
                    for ivarName = 1:size(meanIntTable,2)/2
                        variableNucNamesMeanIntensity = cat(2,variableNucNamesMeanIntensity,{[char(channelNames(ivarName+1,3)) '_' char(channelNames(ivarName+1,4)) '_NucIntensity']});
                        variableCytoNamesMeanIntensity = cat(2,variableCytoNamesMeanIntensity,{[char(channelNames(ivarName+1,3)) '_' char(channelNames(ivarName+1,4)) '_CytoIntensity']});
                        variableNucNamesMedianIntensity = cat(2,variableNucNamesMedianIntensity,{[char(channelNames(ivarName+1,3)) '_' char(channelNames(ivarName+1,4)) '_NucIntensity' ]});
                        variableCytoNamesMedianIntensity = cat(2,variableCytoNamesMedianIntensity,{[char(channelNames(ivarName+1,3)) '_' char(channelNames(ivarName+1,4)) '_CytoIntensity' ]});
                    end
                    variableCytoNamesMeanIntensity = regexprep(variableCytoNamesMeanIntensity, '\s+', '');
                    variableNucNamesMeanIntensity = regexprep(variableNucNamesMeanIntensity, '\s+', '');
                    variableCytoNamesMedianIntensity = regexprep(variableCytoNamesMedianIntensity, '\s+', '');
                    variableNucNamesMedianIntensity = regexprep(variableNucNamesMedianIntensity, '\s+', '');
                else
                    for ivarName = 1:size(meanIntTable,2)/2
                        variableNucNamesMeanIntensity = cat(2,variableNucNamesMeanIntensity,{['Chan_' int2str(ivarName) '_NucIntensity' ]});
                        variableCytoNamesMeanIntensity = cat(2,variableCytoNamesMeanIntensity,{['Chan_' int2str(ivarName) '_CytoIntensity']});
                        variableNucNamesMedianIntensity = cat(2,variableNucNamesMedianIntensity,{['Chan_' int2str(ivarName) '_NucIntensity']});
                        variableCytoNamesMedianIntensity = cat(2,variableCytoNamesMedianIntensity,{['Chan_' int2str(ivarName) '_CytoIntensity']});
                    end
                end
                writetable(array2table([meanIntTable areaTable centroidCellTable],...
                    'VariableNames',[variableNucNamesMeanIntensity variableCytoNamesMeanIntensity...
                    'NucleusArea' 'CytoplasmArea' 'CellPosition_X' 'CellPosition_Y']),...
                    [analysisPath filesep name '_meanNucleiCytoMasked.txt'],'Delimiter','\t')
                
                if isequal(p.MedianIntensity,'true')
                writetable(array2table([medianIntTable areaTable centroidCellTable],...
                    'VariableNames',[variableNucNamesMedianIntensity variableCytoNamesMedianIntensity...
                    'NucleusArea' 'CytoplasmArea' 'CellPosition_X' 'CellPosition_Y']),...
                    [analysisPath filesep name '_medianNucleiCytoMasked.txt'],'Delimiter','\t')
                end
                
            end