function mapAnalysis_Revision(paths,doPlot)

    if nargin < 2 || isempty(doPlot)
        doPlot = false;
    end

    clc
    close all
    drawnow
    fprintf('\n')
    %%% Reliability constaint
    
    warning off all
    if isempty(gcp)
        parpool('local',7);
    end
    pctRunOnAll warning off all
    
    %% Split by animal
    piece = [];
    ag = [];
    spiece = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
        spiece = [spiece; {paths{i}(ind+2:end-4)}];
    end
    upiece = unique(piece);
    
    labels = [{'Saline'} {'CNO'} {'No Injection'}];
    cumDistBins = [0:0.25:5];
    cumDelayDistance = [];
    allDelaySampling = [];
    allDelayTiming = [];
    allNormalTiming = [];
    doExitCoding = false;
    doTwoBack = false;
    doPrediction = false;
    doMedSplit = false;
    doDelay = false;
    doLoc = true;
    doBMRs = false;
    doHallways = false;
    doBehavioralDiagnostics = false;
    velThresh = -2;
    pThresh = 0.05; % 0.05;
    snrThresh = 8;
    delayBins = [0:0.25:5];
    percentPlaceCells = repmat({[]},[1 length(labels)]);
    afw = repmat({[]},[length(labels) length(upiece)]);
    amfr = repmat({[]},[length(labels) length(upiece)]);
    amfr = repmat({[]},[length(labels) length(upiece)]);
    icorr = repmat({[]},[length(labels) length(upiece)]);
    apfshift = repmat({[]},[length(labels) length(upiece)]);
    cellCounts = repmat({[]},[length(labels) length(upiece)]);
    crossAllComp = repmat({[]},[length(labels) length(upiece)]);
    crossExitComp = repmat({[]},[length(labels) length(upiece)]);
    crossTwoBackComp = repmat({[]},[length(labels) length(upiece)]);
    crossAllCompHallway = repmat({[]},[length(labels) length(upiece)]);
    crossAllCompDelayed = repmat({[]},[length(labels) length(upiece)]);
    crossAllMedSplit = repmat({[]},[length(labels) length(upiece)]);
    crossLocComp = repmat({[]},[length(labels) length(upiece)]);
    allPaths = repmat({[]},[length(labels) length(upiece)]);
    allUnmVelBias = repmat({[]},[length(labels) length(upiece)]);
    allMVelBias = repmat({[]},[length(labels) length(upiece)]);
    abmrs = repmat({[]},[length(labels) length(upiece)]);
    afc = repmat({[]},[length(labels) length(upiece)]);
    tmfr = repmat({[]},[length(labels) length(upiece)]);
    tpfr = repmat({[]},[length(labels) length(upiece)]);
    tsi = repmat({[]},[length(labels) length(upiece)]);
    tsnr = repmat({[]},[length(labels) length(upiece)]);
    tsip = repmat({[]},[length(labels) length(upiece)]);
    tcorr = repmat({[]},[length(labels) length(upiece)]);
    tcorrp = repmat({[]},[length(labels) length(upiece)]);
    tcorrpval = repmat({[]},[length(labels) length(upiece)]);
    tcorrpvalp = repmat({[]},[length(labels) length(upiece)]);
    aggVals = repmat({[]},[length(labels) length(upiece)]);
    for mi = 1:length(upiece)
        fprintf(['\n\n\tMouse:  ' num2str(upiece{mi}) '\n']) 
        isM = find(ismember(piece,upiece(mi)));
        for si = 1:length(isM);
            fprintf(['\n\t\tSession:  ' paths{isM(si)}])
            s = load(paths{isM(si)},'processed','properties');
%             s.processed.p = imfilter(s.processed.p,fspecial('gauss',[1 15],2),'same','replicate');
            %id group
            group = find(ismember(labels,paths{isM(si)}(find(ismember(paths{isM(si)},'_'),1,'last')+1:end-4)));
            
            if isempty(group)
                group = 3;
            end
            
            %compute velocity
            vel = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
            vel = imfilter(vel,fspecial('gauss',[1 30],10),'same','replicate');

            cellType = nan(1,4);
            
            %choose cells
            if isfield(s.processed,'exclude')
                inds = (s.processed.splithalf.roomXdoors.p<=pThresh & s.processed.exclude.SFPs);
%                 inds = (s.processed.snr.whole>snrThresh & s.processed.exclude.SFPs);
                
                percentPlaceCells{group} = [percentPlaceCells{group}; nansum(inds)./nansum(s.processed.exclude.SFPs)];
                inc = s.processed.exclude.SFPs;
                cellType = [nansum(inc & s.processed.splithalf.roomXdoors.p>pThresh & s.processed.splithalf.hallwayXdoors.p>pThresh) ...
                    nansum(inc & s.processed.splithalf.roomXdoors.p<=pThresh & s.processed.splithalf.hallwayXdoors.p>pThresh) ...
                    nansum(inc & s.processed.splithalf.roomXdoors.p>pThresh & s.processed.splithalf.hallwayXdoors.p<=pThresh) ...
                    nansum(inc & s.processed.splithalf.roomXdoors.p<=pThresh & s.processed.splithalf.hallwayXdoors.p<=pThresh)];
            else
                inds = (s.processed.splithalf.roomXdoors.p<=pThresh);% & ~s.processed.exclude.SFPs);
%                 inds = s.processed.snr.whole>snrThresh;
                
                percentPlaceCells{group} = [percentPlaceCells{group}; nanmean(inds)];
                if isfield(s.processed.splithalf,'hallwayXdoors')
                    cellType = [nansum(s.processed.splithalf.roomXdoors.p>pThresh & s.processed.splithalf.hallwayXdoors.p>pThresh) ...
                        nansum(s.processed.splithalf.roomXdoors.p<=pThresh & s.processed.splithalf.hallwayXdoors.p>pThresh) ...
                        nansum(s.processed.splithalf.roomXdoors.p>pThresh & s.processed.splithalf.hallwayXdoors.p<=pThresh) ...
                        nansum(s.processed.splithalf.roomXdoors.p<=pThresh & s.processed.splithalf.hallwayXdoors.p<=pThresh)];
                else
                    cellType = [nansum(s.processed.splithalf.roomXdoors.p>pThresh) ...
                        nansum(s.processed.splithalf.roomXdoors.p<=pThresh) 0 0];
                end
            end
            
%             snr = [nanmean(s.calcium.FiltTraces)./nanstd(s.calcium.FiltTraces)];
            
            aggVals(group,mi) = {[aggVals{group,mi}; ...
                [s.processed.snr.whole(inds) s.processed.snr.splithalf(inds,:) ...
                s.processed.splithalf.roomXdoors.p(inds) s.processed.splithalf.wholemap_si.p(inds)]]};
            tsnr(group,mi) = {[tsnr{group,mi}; [s.processed.snr.whole(inds) s.processed.snr.splithalf(inds,:)]]};            
            tsi(group,mi) = {[tsi{group,mi}; s.processed.splithalf.wholemap_si.p]};
            tsip(group,mi) = {[tsip{group,mi}; s.processed.splithalf.wholemap_si.p(inds)]};
            tmfr(group,mi) = {[tmfr{group,mi}; nanmean(s.processed.trace(inds,:),2)]};
            m = mkTraceMaps(s.processed.p,s.processed.trace(inds,:));
            tpfr(group,mi) = {[tpfr{group,mi}; permute(nanmax(nanmax(m,[],1),[],2).*30,[3 1 2])]};
            tcorrp(group,mi) = {[tcorrp{group,mi}; permute(s.processed.splithalf.wholemap.val(1,2,inds),[3 1 2])]};
            tcorr(group,mi) = {[tcorr{group,mi}; permute(s.processed.splithalf.wholemap.val(1,2,:),[3 1 2])]};
            tcorrpvalp(group,mi) = {[tcorrpvalp{group,mi}; s.processed.splithalf.wholemap.p(inds)]};
            tcorrpval(group,mi) = {[tcorrpval{group,mi}; s.processed.splithalf.wholemap.p]};
        
            
            cellCounts{group,mi,1} = [cellCounts{group,mi,1}; {cellType(1)} {cellType(2:4)}];
            str = sprintf(['\n\t\t\tCell Breakdown:  %0.2f%%  %0.2f%%  %0.2f%%  %0.2f%%  '],100.*cellType./nansum(cellType));
            fprintf('%s',str);
            
            gT = s.processed.trace(inds,:);
            
            [isIn isMostRecent] = isInROI(s.processed.p,s.processed.roi.door);
            [isInRoom blah blah2 blah3 indexSinceIn distanceSinceIn] = ...
                isInROI(s.processed.p,s.processed.roi.room);
            indexSinceIn = indexSinceIn.*(1./30);
            half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;
            
            tmp = nan(1,length(cumDistBins));
            tmp(1) = 0;
            for i = 1:length(cumDistBins)-1
                tmp(i+1) = nanmean(distanceSinceIn(indexSinceIn>=cumDistBins(i) & indexSinceIn<cumDistBins(i+1)));
            end
            cumDelayDistance = [cumDelayDistance; tmp];

            allMasks = repmat({[]},[1 4]);
            for i = 1:2
                allMasks{i} = [isMostRecent(i,isInRoom) & half(1,isInRoom) & vel(1,isInRoom)>velThresh];
                allMasks{i+2} = [isMostRecent(i,isInRoom) & ~half(1,isInRoom) & vel(1,isInRoom)>velThresh];
    %             h2 = cumsum(isMostRecent(i,isInRoom));
    %             allMasks{i} = [isMostRecent(i,isInRoom) & h2<h2(end)./2];
    %             allMasks{i+2} = [isMostRecent(i,isInRoom) & ~(h2<h2(end)./2)];
            end
            
            queryMask = false(4);
            queryMask(1:2,3:4) = true;
          
            if ~doLoc
                [map samp allComp] = getMatchedMapsNMasks(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks,queryMask);
            else
                [map samp allComp ivals locvals locivals] = getMatchedMapsNMasks(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks,queryMask);
                crossLocComp{group,mi}  = cat(3,crossLocComp{group,mi},...
                    nanmean(cat(3,locvals(:,:,1,3),locvals(:,:,2,4)),3) - ...
                    nanmean(cat(3,locvals(:,:,2,3),locvals(:,:,1,4)),3));
            end
            allNormalTiming = [allNormalTiming; nansum(samp(:))./30];
            crossAllComp{group,mi} = cat(1,crossAllComp{group,mi},...
                [help_getMaskedVals(allComp(1:2,3:4,:),[true false; false true]) ...
                help_getMaskedVals(allComp(1:2,3:4,:),[false true; true false])]);
            
            [blah1 ivals mfr pfr pfshift fieldCenters] = getMatchedSamplingValues(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks);
            tmp = [help_getMaskedVals(mfr(1:2,3:4,:),[true false; false true]) ...
                help_getMaskedVals(mfr(1:2,3:4,:),[false true; true false])];
            amfr(group,mi) = {[amfr{group,mi}; tmp]};
            icorr(group,mi) = {[icorr{group,mi}; help_getMaskedVals(ivals(1:2,3:4,:),[true false; false true]) ...
                help_getMaskedVals(ivals(1:2,3:4,:),[false true; true false])]};
            
            apfshift(group,mi) = {[apfshift{group,mi}; help_getMaskedVals(pfshift(1:2,3:4,:),[true false; false true]) ...
                help_getMaskedVals(pfshift(1:2,3:4,:),[false true; true false])]};
            
            if doPlot
%                 map = getMatchedMapsNMasks(s.processed.p(:,isInRoom),gT(:,isInRoom),...
%                     [{isMostRecent(1,isInRoom) & vel(1,isInRoom)>velThresh} ...
%                     {isMostRecent(2,isInRoom) & vel(1,isInRoom)>velThresh}]);

                [a b] = sort([tmp(:,2)-tmp(:,1)],'descend');
                
%                 b = b(round([1:length(b)./20:length(b)./2 + length(b)./20]));
                help_plotMaps(map(:,:,b,:),paths(isM(si)));
            end
            
            if doBehavioralDiagnostics
                [unmatchedBehav matchedBehav] = getMatchedMapsBehavioralBiases( ...
                    s.processed.p(:,isInRoom),[{isMostRecent(1,isInRoom) & vel(1,isInRoom)>velThresh} ...
                    {isMostRecent(2,isInRoom) & vel(1,isInRoom)>velThresh}]);
                allUnmVelBias{group,mi} = cat(3,allUnmVelBias{group,mi},...
                    permute(unmatchedBehav(:,:,1,1)-unmatchedBehav(:,:,1,2),[1 2 4 3]));
                allMVelBias{group,mi} = cat(3,allMVelBias{group,mi},...
                    permute(matchedBehav(:,:,1,1)-matchedBehav(:,:,1,2),[1 2 4 3]));
            end
            
            if doMedSplit
                tic
                fprintf('\n\t\t\tComputing Median Split Analysis...  ')
                allDelayMasks = repmat({[]},[1 8]);
                queryMask = false(numel(allDelayMasks));
                queryMask(1:2,3:4) = true;
                queryMask(5:6,7:8) = true;
                for i = 1:2
                    allDelayMasks{i} = [isMostRecent(i,isInRoom) & half(1,isInRoom) & ...
                        vel(1,isInRoom)>velThresh & indexSinceIn(1,isInRoom) < nanmedian(indexSinceIn(isInRoom))];
                    allDelayMasks{i+2} = [isMostRecent(i,isInRoom) & ~half(1,isInRoom) & ...
                        vel(1,isInRoom)>velThresh & indexSinceIn(1,isInRoom) < nanmedian(indexSinceIn(isInRoom))];
                    allDelayMasks{i+4} = [isMostRecent(i,isInRoom) & half(1,isInRoom) & ...
                        vel(1,isInRoom)>velThresh & indexSinceIn(1,isInRoom) >= nanmedian(indexSinceIn(isInRoom))];
                    allDelayMasks{i+6} = [isMostRecent(i,isInRoom) & ~half(1,isInRoom) & ...
                        vel(1,isInRoom)>velThresh & indexSinceIn(1,isInRoom) >= nanmedian(indexSinceIn(isInRoom))];
                end
                
                [map samp allComp] = getMatchedMapsNMasks(s.processed.p(:,isInRoom), ...
                    gT(:,isInRoom),allDelayMasks,queryMask);

                if all(isnan(allComp(:)))
                    continue
                end
                
                delayVals = nan(2,2);
                for q = 1:2
                    delayVals(:,q) = [help_getMaskedVals(allComp(4.*(q-1)+[1:2],4.*(q-1)+[3:4],:),[true false; false true]) ...
                        help_getMaskedVals(allComp(4.*(q-1)+[1:2],4.*(q-1)+[3:4],:),[false true; true false])];
                end

                crossAllMedSplit{group,mi} = cat(3,crossAllMedSplit{group,mi},delayVals);
                
                durat = toc;
                fprintf([num2str(durat) ' sec']);
            end
            
            if doDelay
                tic
                fprintf('\n\t\t\tComputing Delay Analysis...  ')
                allDelayMasks = repmat({[]},[4 length(delayBins)]);
                queryMask = false(numel(allDelayMasks));
                for q = 1:length(delayBins)
                    for i = 1:2
                        allDelayMasks{i,q} = [isMostRecent(i,isInRoom) & half(1,isInRoom) & ...
                            vel(1,isInRoom)>velThresh & indexSinceIn(1,isInRoom) >= delayBins(q)];
                        allDelayMasks{i+2,q} = [isMostRecent(i,isInRoom) & ~half(1,isInRoom) & ...
                            vel(1,isInRoom)>velThresh & indexSinceIn(1,isInRoom) >= delayBins(q)];
                    end
                    queryMask(4.*(q-1)+[1:2],4.*(q-1)+[3:4]) = true;
                end
                os = nan(15,15,length(delayBins));
                for q = 1:length(delayBins)
                    [blah samp] = getMatchedMapsNMasks(s.processed.p(:,isInRoom), ...
                        gT(:,isInRoom),allDelayMasks(:,q)',false(4));
                    os(:,:,q) = samp;
                end
                allDelaySampling = cat(4,allDelaySampling,os);
     
                
                allDelayMasks = allDelayMasks(:)';
                [map samp allComp] = getMatchedMapsNMasks(s.processed.p(:,isInRoom), ...
                    gT(:,isInRoom),allDelayMasks,queryMask);

                if all(isnan(allComp(:)))
                    continue
                end
%                 nansum(samp(:))
                delayVals = nan(2,length(delayBins));
                for q = 1:length(delayBins)
                    delayVals(:,q) = [help_getMaskedVals(allComp(4.*(q-1)+[1:2],4.*(q-1)+[3:4],:),[true false; false true]) ...
                        help_getMaskedVals(allComp(4.*(q-1)+[1:2],4.*(q-1)+[3:4],:),[false true; true false])];
                end

                allDelayTiming = [allDelayTiming; permute(nansum(nansum(os,1),2),[1 3 2])./30  ];
                
%                 if nansum(samp(:))./30 < 30 % Must include at least 30s of data per comparison to be included
%                     continue
%                 end
                
                plot(indexSinceIn(:,isInRoom))
                crossAllCompDelayed{group,mi} = cat(3,crossAllCompDelayed{group,mi},delayVals);

                durat = toc;
                fprintf([num2str(durat) ' sec']);

            end
            
            allPaths{group,mi} = [allPaths{group,mi} [{s.processed.p(:,~isInRoom)}; ...
                {s.processed.p(:,isInRoom&isMostRecent(1,:))}; {s.processed.p(:,isInRoom&isMostRecent(2,:))}]];
            if doBMRs
                [blah1 blah2 blah3 blah4 blah5 fieldCenters bmrs] = ...
                    getMatchedSamplingValues(s.processed.p(:,isInRoom),gT(:,isInRoom), ...
                    [{isMostRecent(1,isInRoom)} {isMostRecent(2,isInRoom)}]);
                abmrs{group,mi} = [abmrs{group,mi}; permute(bmrs(1,2,:),[1 3 2])];
                afc{group,mi} = [afc{group,mi}; fieldCenters(:,:,1) fieldCenters(:,:,2)];
            end
           
            if doTwoBack
                [isIn isMostRecent] = isInROI(s.processed.p,(s.processed.roi.door));
                isMostRecent(:,~isInRoom) = false;
                isSecondMostRecent = nthMostRecent(isMostRecent,3);
                allMasks = repmat({[]},[1 4]);
                for i = 1:2
                    allMasks{i} = [isSecondMostRecent(i,isInRoom) & half(1,isInRoom) & vel(1,isInRoom)>velThresh];
                    allMasks{i+2} = [isSecondMostRecent(i,isInRoom) & ~half(1,isInRoom) & vel(1,isInRoom)>velThresh];
                end
                queryMask = false(4);
                queryMask(1:2,3:4) = true;

                [map samp twoBackComp] = getMatchedMapsNMasks( ...
                    s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks,queryMask);
                crossTwoBackComp{group,mi} = [crossTwoBackComp{group,mi}; ...
                    [help_getMaskedVals(twoBackComp(1:2,3:4,:),[true false; false true]) ...
                    help_getMaskedVals(twoBackComp(1:2,3:4,:),[false true; true false])]];
                
            end
            
            if doExitCoding
                [isIn isMostRecent] = isInROI(fliplr(s.processed.p),(s.processed.roi.door));
                isMostRecent = fliplr(isMostRecent);
                allMasks = repmat({[]},[1 4]);
                for i = 1:2
                    allMasks{i} = [isMostRecent(i,isInRoom) & half(1,isInRoom) & vel(1,isInRoom)>velThresh];
                    allMasks{i+2} = [isMostRecent(i,isInRoom) & ~half(1,isInRoom) & vel(1,isInRoom)>velThresh];
                end

                queryMask = false(4);
                queryMask(1:2,3:4) = true;

                [map samp exitComp] = getMatchedMapsNMasks( ...
                    s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks,queryMask);
                crossExitComp{group,mi} = [crossExitComp{group,mi}; ...
                    [help_getMaskedVals(exitComp(1:2,3:4,:),[true false; false true]) ...
                    help_getMaskedVals(exitComp(1:2,3:4,:),[false true; true false])]];
            end

            if doHallways
                if isfield(s.processed,'exclude')
                    hinds = (s.processed.splithalf.hallwayXdoors.p<=pThresh & s.processed.exclude.SFPs);
                else
                    hinds = (s.processed.splithalf.hallwayXdoors.p<=pThresh);
                end
                hgT = s.processed.trace(hinds,:);
                isIn = inpolygon(s.processed.p(1,:)',s.processed.p(2,:)',...
                    s.processed.roi.hallway(:,1),s.processed.roi.hallway(:,2));        
                cp = help_collapseToLine(s.processed.p,s.processed.roi.hallway_linear);

                isGood = isIn; %%% Only include full movements
                isDir = false(2,length(isIn));
                oneEntrance = cp > [nanmax(cp(isIn))-nanmin(cp(isIn))]./2;
                while any(isGood)
                    start = find(isGood,1,'first');
                    stop = find(~isGood(start:end),1,'first')-2;
                    if isempty(stop)
                        stop = length(isGood)-start;
                    end
        %             if oneEntrance(start)~=oneEntrance(start+stop)
                        if oneEntrance(start)
                            isDir(1,start:start+stop) = true;
                        else
                            isDir(2,start:start+stop) = true;
                        end
        %             end
                    isGood(start:start+stop) = false;
                end
                
                [maps samp allComp ival] = getMatchedLinMapsNMasks(cp(isIn),hgT(:,isIn),...
                    [{isDir(1,isIn) & half(isIn)} {isDir(2,isIn)& half(isIn)} ...
                    {isDir(1,isIn) & ~half(isIn)} {isDir(2,isIn)& ~half(isIn)}]);

                crossAllCompHallway{group,mi} = cat(1,crossAllCompHallway{group,mi},...
                    [help_getMaskedVals(allComp(1:2,3:4,:),[true false; false true]) ...
                    help_getMaskedVals(allComp(1:2,3:4,:),[false true; true false])]);
                

            end
        end
    end
    close all
    
    eliminate = all(cellfun(@isempty,crossAllComp),2);
    labels(eliminate) = [];
    amfr(eliminate,:) = [];
    icorr(eliminate,:) = [];
    afw(eliminate,:) = [];
    cellCounts(eliminate,:) = [];
    crossAllComp(eliminate,:) = [];
    crossAllCompDelayed(eliminate,:) = [];
    allPaths(eliminate,:) = [];
    crossLocComp(eliminate,:) = [];
    allUnmVelBias(eliminate,:) = [];
    allMVelBias(eliminate,:) = [];
    apfshift(eliminate,:) = [];
    crossAllMedSplit(eliminate,:) = [];
    abmrs(eliminate,:) = [];
    afc(eliminate,:) = [];
    tsnr(eliminate,:) = [];
    tmfr(eliminate,:) = [];
    tpfr(eliminate,:) = [];
    tsi(eliminate,:) = [];
    tsip(eliminate,:) = [];
    tcorr(eliminate,:) = [];
    tcorrp(eliminate,:) = [];
    tcorrpval(eliminate,:) = [];
    tcorrpvalp(eliminate,:) = [];
    crossAllCompHallway(eliminate,:) = [];
    crossExitComp(eliminate,:) = [];
    crossTwoBackComp(eliminate,:) = [];
    aggVals(eliminate,:) = [];
    
    
    slashInds = find(ismember(paths{1},'/'));
    root = ['Plots/REVISION/Summary' paths{1}(slashInds(1):slashInds(2)-1)];
    
    if doLoc
        % Effect size as a function of distance to nearest entrance
        
        tmp = nanmedian(cat(3,crossLocComp{:}),3);
        entrances = [7 13; 1 6];
        [x y] = meshgrid(1:15,1:15);
        d2e1 = sqrt([[x-entrances(1,2)].^2 + [y-entrances(1,1)].^2]);
        d2e2 = sqrt([[x-entrances(2,2)].^2 + [y-entrances(2,1)].^2]);
        d2ne = nanmin(d2e1,d2e2);
        rd2ne = repmat(d2ne,[1 1 length(tmp(1,1,:))]);
        a = rd2ne(:);
        b = tmp(:);
        a(isnan(b)) = [];
        b(isnan(b)) = [];
        [r pval] = corr(a,b)
        
        tmp = cat(3,crossLocComp{:});
        entrances = [7 13; 1 6];
        [x y] = meshgrid(1:15,1:15);
        d2e1 = sqrt([[x-entrances(1,2)].^2 + [y-entrances(1,1)].^2]);
        d2e2 = sqrt([[x-entrances(2,2)].^2 + [y-entrances(2,1)].^2]);
        d2ne = nanmin(d2e1,d2e2);
        rd2ne = repmat(d2ne,[1 1 length(tmp(1,1,:))]);
        a = rd2ne(:);
        b = tmp(:);
        a(isnan(b)) = [];
        b(isnan(b)) = [];
        [r pval] = corr(a,b)
    end
    
    
%     figure
%     set(gcf,'position',[50 50 300 300])
%     subplot(2,1,1)
%     hist(amfr{1}(:,2)-amfr{1}(:,1),[-0.75:0.05:0.75]);
%     hold on
%     plot([nanmedian(amfr{1}(:,2)-amfr{1}(:,1)) nanmedian(amfr{1}(:,2)-amfr{1}(:,1))],get(gca,'ylim'),...
%         'linestyle','-','color','r','linewidth',2)
%     subplot(2,1,2)
%     hist(amfr{2}(:,2)-amfr{2}(:,1),[-0.75:0.05:0.75]);
%     hold on
%     plot([nanmedian(amfr{2}(:,2)-amfr{2}(:,1)) nanmedian(amfr{2}(:,2)-amfr{2}(:,1))],get(gca,'ylim'),...
%         'linestyle','-','color','r','linewidth',2)
%     saveFig(gcf,[root '/RR_Example'],[{'pdf'} {'tiff'}]);
    
% % % %     tmp = nanmean(allDelaySampling,4);
% % % %     figure
% % % %     set(gcf,'position',[50 50 700 500])
% % % %     for i = 1:2
% % % %         subplot(1,2,i)
% % % %         imagesc(tmp(1:13,1:13,i)./30)
% % % %         colormap parula
% % % %         alpha(double(0~=(tmp(1:13,1:13,i)./30)))
% % % %         caxis([0 2])
% % % %         colorbar
% % % %         axis equal
% % % %         axis off
% % % %     end
% % % %     saveFig(gcf,[root '/Sampling_With_Delay'],[{'pdf'} {'tiff'}])

    
%     mkWhisker(cumDelayDistance,cumDistBins)

% % %     step = 1;
% % %     for ui = 1:length(upiece)
% % %         for i = 1:length(labels)
% % %             figure
% % %             set(gcf,'position',[500 450 300.*length(allPaths{i,ui}) length(labels).*300])
% % %             for j = 1:length(allPaths{i,ui}(1,:))
% % %                 subplot(1,length(allPaths{i,ui}),j)
% % % %                 s1 = scatter(allPaths{i,ui}{1,j}(1,1:step:end),allPaths{i,ui}{1,j}(2,1:step:end),3,'k',...
% % % %                     'markerfacecolor','k');
% % % %                 hold on
% % % %                 s2 = scatter(allPaths{i,ui}{2,j}(1,1:step:end),allPaths{i,ui}{2,j}(2,1:step:end),3,[0.3 0.3 0.9],...
% % % %                     'markerfacecolor',[0.3 0.3 0.9]);
% % % %                 s3 = scatter(allPaths{i,ui}{3,j}(1,1:step:end),allPaths{i,ui}{3,j}(2,1:step:end),3,[0.3 0.8 0.9],...
% % % %                     'markerfacecolor',[0.3 0.8 0.9]);
% % % %                 plot(allPaths{i,ui}{1,j}(1,1:step:end),allPaths{i,ui}{1,j}(2,1:step:end),'linestyle','none',...
% % % %                     'markerfacecolor','k','marker','o','markersize',2.5,'markeredgecolor','none');
% % % %                 hold on
% % % %                 plot(allPaths{i,ui}{2,j}(1,1:step:end),allPaths{i,ui}{2,j}(2,1:step:end),'linestyle','none',...
% % % %                     'markerfacecolor',[0.3 0.3 0.9],'marker','o','markersize',2.5,'markeredgecolor','none');
% % % %                 plot(allPaths{i,ui}{3,j}(1,1:step:end),allPaths{i,ui}{3,j}(2,1:step:end),'linestyle','none',...
% % % %                     'markerfacecolor',[0.3 0.8 0.9],'marker','o','markersize',2.5,'markeredgecolor','none');
% % %                 
% % %                 breaks = [1 find(nansum(diff(allPaths{i,ui}{1,j},[],2).^2)>100) length(allPaths{i,ui}{1,j})];
% % %                 hold on
% % %                 for bi = 1:length(breaks)-1
% % %                     plot(allPaths{i,ui}{1,j}(1,breaks(bi)+1:breaks(bi+1)),...
% % %                         allPaths{i,ui}{1,j}(2,breaks(bi)+1:breaks(bi+1)),'linestyle','-',...
% % %                         'color','k');
% % %                 end
% % %                 breaks = [1 find(nansum(diff(allPaths{i,ui}{2,j},[],2).^2)>100) length(allPaths{i,ui}{2,j})];
% % %                 for bi = 1:length(breaks)-1
% % %                     plot(allPaths{i,ui}{2,j}(1,breaks(bi)+1:breaks(bi+1)),...
% % %                         allPaths{i,ui}{2,j}(2,breaks(bi)+1:breaks(bi+1)),'linestyle','-',...
% % %                         'color',[0.3 0.3 0.9]);
% % %                 end
% % %                 breaks = [1 find(nansum(diff(allPaths{i,ui}{3,j},[],2).^2)>100) length(allPaths{i,ui}{3,j})];
% % %                 for bi = 1:length(breaks)-1
% % %                     plot(allPaths{i,ui}{3,j}(1,breaks(bi)+1:breaks(bi+1)),...
% % %                         allPaths{i,ui}{3,j}(2,breaks(bi)+1:breaks(bi+1)),'linestyle','-',...
% % %                         'color',[0.3 0.8 0.9]);
% % %                 end
% % % %                 hold on
% % % %                 plot(allPaths{i,ui}{2,j}(1,1:step:end),allPaths{i,ui}{2,j}(2,1:step:end),'linestyle','-',...
% % % %                     'color',[0.3 0.3 0.9]);
% % % %                 plot(allPaths{i,ui}{3,j}(1,1:step:end),allPaths{i,ui}{3,j}(2,1:step:end),'linestyle','-',...
% % % %                     'color',[0.3 0.8 0.9]);
% % % 
% % %                 axis equal
% % %                 axis off
% % %             end
% % %             tmp = upiece{ui};
% % %             slashInds = find(ismember(tmp,'/'));
% % %             saveFig(gcf,[root '/Paths/' upiece{ui}(slashInds(end)+1:end) '_' labels{i}],[{'pdf'} {'tiff'}]);
% % %             close(gcf)
% % %         end
% % %     end

    
    
    
    h1 = mkPie(cellCounts);
%     lgd = legend([h1 h2],[{'Nonspatial'} {'Spatial'} {'Room'} {'Hallway'} {'Both'}],...
%         'box','off','location','northoutside','orientation','horizontal');
    saveFig(gcf,[root '/CellClassification'],[{'pdf'} {'tiff'}]);
    
    figure
    set(gcf,'position',[50 50 200.*length(tsnr(1,:)).*length(tsnr(:,1)) 300.*1])
    for mi = 1:length(tsnr(1,:))
        subplot(1,length(tsnr(1,:)),mi)

        tmp = upiece{mi};
        slashInds = find(ismember(tmp,'/'));

        tmp2 = tsnr(:,mi);
        doPlot = repmat({[]},[length(tmp2(:,1)) 3]);
        for k = 1:length(tmp2(:,1))
            for hi = 1:3
                doPlot{k,hi} = tmp2{k}(tmp2{k}(:,hi)~=inf,hi);
            end
        end
        doPlot = doPlot';
%         mkGraph(doPlot',labels, ...
%             [{' Whole Trial'} {' 1st Half'} {' 2nd Half'}]);
        mkWhisker([doPlot(:)]',[{' Whole Trial'} {' 1st Half'} {' 2nd Half'}])
        set(gca,'xticklabelrotation',45)
        ylabel('SNR')
%         set(gca,'ylim',[0 15]);
        title(upiece{mi}(slashInds(end)+1:end))
    end
    saveFig(gcf,[root '/SNR'],[{'pdf'} {'tiff'}]);
    
    
    figure
    set(gcf,'position',[50 450 175.*length(upiece) 200])
    toPlot = repmat({[]},[2 length(labels) length(upiece)]);
    lim = cat(1,crossAllComp{:});
    lim = ceil(nanmax(lim(:)).*10)./10;
    for ui = 1:length(upiece)
        for group = 1:length(labels)
            toPlot{1,group,ui} = crossAllComp{group,ui}(:,1);
            toPlot{2,group,ui} = crossAllComp{group,ui}(:,2);
        end
        subplot(1,length(upiece),ui)
        mkGraph(toPlot(:,:,ui),{'Entryway'},[{'Same'} {'Diff.'}])
        ylabel('Correlation (r)')
        tmp = upiece{ui};
        slashInds = find(ismember(tmp,'/'));
        title(upiece{ui}(slashInds(end)+1:end))
        set(gca,'ylim',[0 lim])
    end
    saveFig(gcf,[root '/PopVecRemapping_AnimalWise'],[{'pdf'} {'tiff'}]);
    
    
    vals = nan([nanmax(cellfun(@length,crossAllComp(:))) size(crossAllComp')]);
    for ui = 1:length(upiece)
        for group = 1:length(labels)
            vals(1:length(crossAllComp{group,ui}(:,1)),ui,group) = crossAllComp{group,ui}(:,1)-crossAllComp{group,ui}(:,2);
        end
    end
%     timePlot = repmat({[]},[length(labels) length(vals(:,1,1))]);
%     for i = 1:length(vals(:,1,1))
%         for k = 1:length(vals(1,1,:))
%             timePlot{k,i} = vals(i,:,k)'; 
%         end
%     end
%     fid = figure;
%     set(gcf,'position',[50 450 350 200])
%     mkGraph(timePlot')
%     set(gca,'ylim',[-0.1 0.2])
%     saveFig(gcf,[root '/PopVecRemappingXSession'],[{'pdf'} {'tiff'}]);
%     close(fid);
%     
%     timePlot = [{reshape(vals(1:4,:,1),1,[])'} {reshape(vals(5:9,:,1),1,[])'}; ...
%         {reshape(vals(1:4,:,2),1,[])'} {reshape(vals(5:9,:,2),1,[])'}];
%     fid = figure;
%     set(gcf,'position',[50 450 175 200])
%     mkGraph(timePlot')
%     set(gca,'ylim',[-0.1 0.2])
%     saveFig(gcf,[root '/PopVecRemappingXSession_Binned'],[{'pdf'} {'tiff'}]);
%     close(fid);
%     
%     timePlot = [timePlot(1,:)'; timePlot(2,:)'];
%     stats = nan(length(timePlot));
%     pval = nan(length(timePlot));
%     istats = nan(size(timePlot));
%     ipval = nan(size(timePlot));
%     for i = 1:length(timePlot(:,1))
%         [pv h ss] = signrank(timePlot{i});
%         ipval(i) = pv;
%         for j = i+1:length(timePlot(:,1))
%             [pv h ss] = ranksum(timePlot{i},timePlot{j});
%             pval(i,j) = pv;
% %             stats(i) = ss.zval;
%         end
%     end
    
%     vals = nan(22,4);
%     for ui = 1:length(upiece)
%         vals(1:length(crossAllComp{1,ui}(:,1)),ui) = crossAllComp{1,ui}(:,1)-crossAllComp{1,ui}(:,2);
%     end
%     timePlot = [{vals(1:7,:)} {vals(8:14,:)} {vals(15:end,:)}];
%     for i = 1:3
%         timePlot{i} = timePlot{i}(~isnan(timePlot{i}));
%         [pval h stat] = signrank(timePlot{i})
%     end
%     apv = nan(3);
%     astat = nan(3);
%     for i = 1:3
%         for j = i+1:3
%             [pval h stat] = ranksum(timePlot{i},timePlot{j});
%             apv(i,j) = pval;
%             astat(i,j) = stat.zval;
%         end
%     end
%     mkGraph(timePlot)
    
    figure
    set(gcf,'position',[50 450 175 200])
    toPlot = repmat({[]},[2 length(labels)]);
    diffPlot = repmat({[]},[1 length(labels)]);
    for group = 1:length(labels)
        tmp = cat(1,crossAllComp{group,:});
        toPlot{1,group} = tmp(:,1);
        toPlot{2,group} = tmp(:,2);
        diffPlot{group} = tmp(:,1)-tmp(:,2);
    end
    mostRecent = toPlot;
    aDiff = diffPlot;
    mkGraph(toPlot)
    ylabel('Correlation (r)')
    set(gca,'ylim',[0 1])
    saveFig(gcf,[root '/PopVecRemapping_Combined'],[{'pdf'} {'tiff'}]);
    
    figure
    set(gcf,'position',[50 450 175 200])
    mkGraph(diffPlot)
    ylabel('Correlation (r)')
    set(gca,'ylim',[-0.4 0.4])
    saveFig(gcf,[root '/PopVecRemapping_Combined_Difference'],[{'pdf'} {'tiff'}]);
      
    if doExitCoding
        figure
        set(gcf,'position',[50 450 175 200])
        toPlot = repmat({[]},[2 length(labels)]);
        diffPlot = repmat({[]},[1 length(labels)]);
        for group = 1:length(labels)
            tmp = cat(1,crossExitComp{group,:});
            toPlot{1,group} = tmp(:,1);
            toPlot{2,group} = tmp(:,2);
            diffPlot{group} = tmp(:,1)-tmp(:,2);
        end
        bDiff = diffPlot;
        exitCoding = toPlot;
        mkGraph(toPlot)
        ylabel('Correlation (r)')
        set(gca,'ylim',[0 lim])
        saveFig(gcf,[root '/ExitCoding/PopVecRemapping_Combined'],[{'pdf'} {'tiff'}]);

        figure
        set(gcf,'position',[50 450 175 200])
        mkGraph(diffPlot)
        ylabel('Correlation (r)')
    %     set(gca,'ylim',[0 lim])
        saveFig(gcf,[root '/ExitCoding/PopVecRemapping_Combined_Difference'],[{'pdf'} {'tiff'}]);
    end
    
    if doTwoBack
        figure
        set(gcf,'position',[50 450 175 200])
        toPlot = repmat({[]},[2 length(labels)]);
        diffPlot = repmat({[]},[1 length(labels)]);
        for group = 1:length(labels)
            tmp = cat(1,crossTwoBackComp{group,:});
            toPlot{1,group} = tmp(:,1);
            toPlot{2,group} = tmp(:,2);
            diffPlot{group} = tmp(:,1)-tmp(:,2);
        end
        twoBack = toPlot;
        cDiff = diffPlot;
        mkGraph(toPlot)
        ylabel('Correlation (r)')
        set(gca,'ylim',[0 lim])
        saveFig(gcf,[root '/TwoBackCoding/PopVecRemapping_Combined'],[{'pdf'} {'tiff'}]);

        figure
        set(gcf,'position',[50 450 175 200])
        mkGraph(diffPlot)
        ylabel('Correlation (r)')
    %     set(gca,'ylim',[0 lim])
        saveFig(gcf,[root '/TwoBackCoding/PopVecRemapping_Combined_Difference'],[{'pdf'} {'tiff'}]);
    end
    
    if doTwoBack && doExitCoding
        figure
        set(gcf,'position',[50 450 175 200])
        mkGraph([twoBack mostRecent exitCoding],[{'Two Back'} {'Most Recent'} {'Exit Coding'}]);
        saveFig(gcf,[root '/MultipleCodings/PopVecRemapping_Combined'],[{'pdf'} {'tiff'}]);
        figure
        set(gcf,'position',[50 450 175 200])
        mkGraph([cDiff aDiff bDiff],[{'Two Back'} {'Most Recent'} {'Exit Coding'}]);
        saveFig(gcf,[root '/MultipleCodings/PopVecRemapping_Combined_Difference'],[{'pdf'} {'tiff'}]);
        tmp = [cDiff aDiff bDiff];
        pvals = nan(3);
        astats = nan(3);
        for i = 1:3
            [pval h stats] = signrank(tmp{i})
            for j = i+1:3
                [pval h stats] = ranksum(tmp{i},tmp{j});
                pvals(i,j) = pval;
                astats(i,j) = stats.zval;
            end
        end
    end
    
    if doHallways
        figure
        set(gcf,'position',[50 450 175 200])
        toPlotH = repmat({[]},[2 length(labels)]);
        diffPlotH = repmat({[]},[1 length(labels)]);
        for group = 1:length(labels)
            tmp2 = cat(1,crossAllCompHallway{group,:});
            toPlotH{1,group} = tmp2(:,1);
            toPlotH{2,group} = tmp2(:,2);
            diffPlotH{group} = tmp2(:,1)-tmp2(:,2);
        end
        mkGraph(toPlotH)
        ylabel('Correlation (r)')
        set(gca,'ylim',[0 1])
        saveFig(gcf,[root '/Hallway_PopVecRemapping_Combined'],[{'pdf'} {'tiff'}]);

        figure
        set(gcf,'position',[50 450 175 200])
        mkGraph(diffPlotH)
        ylabel('Correlation (r)')
        set(gca,'ylim',[-0.2 1])
        saveFig(gcf,[root '/Hallway_PopVecRemapping_Combined_Difference'],[{'pdf'} {'tiff'}]);
        
        figure
        set(gcf,'position',[50 50 250.*length(labels) 250])
        for group = 1:length(labels)
            tmp = cat(1,crossAllComp{group,:});
            tmp2 = cat(1,crossAllCompHallway{group,:});
            subplot(1,length(labels),group)
            scatter(tmp(:,1)-tmp(:,2),tmp2(:,1)-tmp2(:,2))
            ylabel('Hallway (same - diff.)')
            xlabel('Room (same - diff.)')
            [r pval] = corr(tmp(:,1)-tmp(:,2),tmp2(:,1)-tmp2(:,2));
            text(nanmin(get(gca,'xlim'))+0.05,nanmax(get(gca,'ylim')+0.05), ...
                sprintf('r = %0.3f, p = %0.3f',[r pval]),'horizontalalignment','Left',...
                                'verticalalignment','middle','fontname','arial',...
                                'fontweight','bold','fontsize',9)
            lsline
            axis square
            set(gca,'ylim',[-0.2 1])
            set(gca,'xlim',[-0.1 0.3])
        end
        saveFig(gcf,[root '/Hallway_vs_Room_PopVecRemapping'],[{'pdf'} {'tiff'}]);
    end
    
    if doBMRs
        for group = 1:length(labels)
            tmp = cat(1,abmrs{group,:});
            figure
            set(gcf,'position',[50 450 175 200])
            mkGraph(tmp)
            ylabel('Correlation (r)')
            set(gca,'ylim',[-0.25 1])
            saveFig(gcf,[root '/BestMatchRotations_' labels{group}],[{'pdf'} {'tiff'}]);
        end
        for group = 1:length(labels)
            tmp = cat(1,afc{group,:});
            [r1 pval1] = corr(tmp(:,1),tmp(:,3));
            [r2 pval2] = corr(tmp(:,2),tmp(:,4));
            figure
            set(gcf,'position',[50 450 400 200])
            subplot(1,2,1)
            scatter(tmp(:,1),tmp(:,3),2)
            xlabel('Entryway A')
            ylabel('Entryway B')
            hold on
            plot([0 15],[0 15],'color','k','linewidth',0.5);
            axis equal
            axis square
            set(gca,'xlim',[0 12],'ylim',[0 12])
            subplot(1,2,2)
            scatter(tmp(:,2),tmp(:,4),2)
            xlabel('Entryway A')
            ylabel('Entryway B')
            hold on
            plot([0 15],[0 15],'color','k','linewidth',0.5);
            axis equal
            axis square
            set(gca,'xlim',[0 12],'ylim',[0 12])
            saveFig(gcf,[root '/FieldCOMCorrelation_' labels{group}],[{'pdf'} {'tiff'}]);
        end
    end
    
    figure
    set(gcf,'position',[500 450 300.*length(upiece) length(labels).*300])
    for ui = 1:length(upiece)
        for i = 1:length(labels)
            subplot(length(labels),length(upiece),(i-1).*length(upiece)+ui)
            h = cumHist(amfr{i,ui},[0:0.025:1]);
            xlabel(sprintf('Split-half change\nin mean firing rate (%%)'),'interpreter','none')
            ylabel('Cumulative proportion')
            legend([h{1}(1) h{2}(1)],[{'Same'} ...
                {'Diff.'} {'Shuffled'}],'location','southeast','fontname','arial',...
                'fontsize',9,'fontweight','bold','color','none','box','off')
            set(gca,'xticklabel',100.*cellfun(@str2num,get(gca,'xticklabel')))
            tmp = upiece{ui};
            slashInds = find(ismember(tmp,'/'));
            title({upiece{ui}(slashInds(end)+1:end)})
            axis square
        end
    end
    drawnow
    saveFig(gcf,[root '/RateRemapping_AnimalWise'],[{'tiff'} {'pdf'}]);
    
    figure
    set(gcf,'position',[500 450 300 length(labels).*300])
    for i = 1:length(labels)
        subplot(length(labels),1,i)
        h = cumHist(cat(1,amfr{i,:}),[0:0.025:1]);
%         h = cumHist(cat(1,amfr{i,:}),[0:0.025:1]);
        xlabel('Mean Rate Change (%)')
        ylabel('Count')
        legend([h{1}(1) h{2}(1)],[{'Within'} ...
            {'Across'} {'Shuffled'}],'location','southeast','fontname','arial',...
            'fontsize',9,'fontweight','bold','color','none','box','off')
        set(gca,'xticklabel',100.*cellfun(@str2num,get(gca,'xticklabel')))
        axis square
    end
    drawnow
    saveFig(gcf,[root '/RateRemapping_Combined'],[{'tiff'} {'pdf'}]);
    
    tmp = cat(1,apfshift{i,:});
    [pval h stat] = signrank(tmp(:,1),tmp(:,2));
    
    figure
    set(gcf,'position',[500 450 300 length(labels).*300])
    for i = 1:length(labels)
        subplot(length(labels),1,i)
        h = cumHist(cat(1,apfshift{i,:}).*2.5,[0:0.025:15.*2.5]);
%         h = cumHist(cat(1,amfr{i,:}),[0:0.025:1]);
%         xlabel('Individual Cell Correlations')
        ylabel('Count')
        legend([h{1}(1) h{2}(1)],[{'Within'} ...
            {'Across'} {'Shuffled'}],'location','southeast','fontname','arial',...
            'fontsize',9,'fontweight','bold','color','none','box','off')
        set(gca,'xticklabel',cellfun(@str2num,get(gca,'xticklabel')))
        axis square
    end
    drawnow
    saveFig(gcf,[root '/LocationRemapping_Combined'],[{'tiff'} {'pdf'}]);

    
    if doBehavioralDiagnostics
        figure
        set(gcf,'position',[500 450 400 length(labels).*300])
        locBiasVals = repmat({[]},[1 length(labels)]);
        for i = 1:length(labels)
            subplot(length(labels),1,i)
            plotVals = cat(3,allMVelBias{i,:});
            mask = nansum(~isnan(plotVals),3)<5;
            plotVals(repmat(mask,[1 1 length(plotVals(1,1,:))])) = nan;
            imagesc(nanmean(plotVals,3))
            alpha(double(~isnan(nanmean(plotVals,3))))
            colormap jet
            caxis([-20 20])
            colorbar
            axis equal
            axis off
            
            tmp = reshape(plotVals,[numel(plotVals(:,:,1)) length(plotVals(1,1,:))]);
            locBiasVals{1} = tmp;
            apv = nan(length(tmp(:,1)),1);
            for j = 1:length(tmp(:,1))
%                 [h pval ci tstat] = ttest(tmp(j,~isnan(tmp(j,:)))',0);
                if all(isnan(tmp(j,:)))
                    continue
                end
                [pval h stat] = signrank(tmp(j,~isnan(tmp(j,:)))');
                apv(j) = pval;
            end
            apv = reshape(apv,size(plotVals(:,:,1)));

            for j = 1:length(plotVals(1,:,1))
                for k = 1:length(plotVals(:,1,1))
                    if apv(j,k) < 0.05
                        text(k,j+0.075,'*','horizontalalignment','center',...
                            'verticalalignment','middle','fontname','arial',...
                            'fontweight','bold','fontsize',9)
                    end
                end
            end
        end
        
        drawnow
        saveFig(gcf,[root '/VelocityBiasXLocation_Matched'],[{'tiff'} {'pdf'}]);
    end
    
    if doMedSplit
        figure
        set(gcf,'position',[50 450 175 200])
        toPlot = repmat({[]},[2 length(labels)]);
        diffPlot = repmat({[]},[1 length(labels)]);
        for group = 1:length(labels)
            tmp = cat(3,crossAllMedSplit{group,:});
            tmp = permute(tmp,[3 2 1]);
            tmp = -diff(tmp,[],3);
            toPlot{1,group} = tmp(:,1);
            toPlot{2,group} = tmp(:,2);
            diffPlot{group} = tmp(:,1)-tmp(:,2);
        end
        mkGraph(toPlot)
        ylabel('Correlation (r)')
        saveFig(gcf,[root '/PopVecRemapping_MedianSplit_Combined_Difference'],[{'pdf'} {'tiff'}]);
        [pval h stat] = signrank(toPlot{1});
        [pval h stat] = signrank(toPlot{2});
        [pval h stat] = ranksum(toPlot{1},toPlot{2});
    end
    
    
    figure
    set(gcf,'position',[500 450 300 length(labels).*200])
    h = mkLine({cumDelayDistance},cumDistBins);
    xlabel('Time since entry (s)')
    ylabel('Cumulative Distance (cm)')
    set(gca,'ylim',[0 80]);
    set(gca,'xtick',[0:1:delayBins(end)],'xticklabel',[0:1:delayBins(end)])
    saveFig(gcf,[root '/DelayedRemapping_CumulativeDistanceTraveled_' ...
        num2str(cumDistBins(1)) '_to_' num2str(cumDistBins(end))],[{'tiff'} {'pdf'}]);
    
    if doDelay
        figure
        set(gcf,'position',[500 450 300 length(labels).*200])
        for i = 1:length(labels)
            tmp = cat(3,crossAllCompDelayed{i,:});
            tmp = permute(tmp,[3 2 1]);
            tmp = tmp(allDelayTiming(:,end)>15,:,:);
            subplot(length(labels),1,i)
            h = mkLine([{tmp(:,:,1)-tmp(:,:,2)}],delayBins);
%             legend(h,[{'Same entry'} {'Other entry'} {'Difference'}],'location','northeast')
    %         h = cumHist(cat(1,amfr{i,:}),[0:0.025:1]);
            xlabel('Time since entry (s)')
            ylabel('Same - different')
%             legend([h(1)],[{'Same - different'}],'location','eastoutside','fontname','arial',...
%                 'fontsize',9,'fontweight','bold','color','none','box','off')
            set(gca,'ylim',[-0.1 0.2])
%             set(gca,'ylim',[round(10.*nanmin(nanmin(get(gca,'ylim')),0))./10 0.1+round(10.*nanmax(get(gca,'ylim')))./10])
            set(gca,'xtick',[0:1:delayBins(end)],'xticklabel',[0:1:delayBins(end)])
%             [h pval ci tstat] = ttest(tmp(:,:,1)-tmp(:,:,2),0);
            for i = 1:length(tmp(1,:,1))
                [pval(i) h stat] = signrank(tmp(:,i,1)-tmp(:,i,2));
            end
            if any(pval<0.05)
                plot(delayBins(pval<0.05),nanmax(nanmax(tmp(:,:,1)-tmp(:,:,2)))+0.02, ...
                    'marker','o','color','k','markerfacecolor','k','markersize',1)
%                 plot(delayBins(pval<0.05),0.075, ...
%                     'marker','o','color','k','markerfacecolor','k','markersize',1)
            end
        end
        
        saveFig(gcf,[root '/DelayedRemapping_Combined_' ...
            num2str(delayBins(1)) '_to_' num2str(delayBins(end))],[{'tiff'} {'pdf'}]);
        
        figure
        set(gcf,'position',[50 50 400 250])
        plot(allDelayTiming','color',[0.75 0.75 0.9])
        hold on
        plot(nanmean(allDelayTiming),'color',[0.1 0.1 0.9],'linewidth',1.5)
        set(gca,'xtick',[1:delayBins(end)],'xticklabel',delayBins(get(gca,'xtick')), ...
            'xlim',[0.75 length(delayBins)+0.25])
        set(gca,'ylim',[0 120])
        saveFig(gcf,[root '/DelayedRemapping_Sampling_' ...
            num2str(delayBins(1)) '_to_' num2str(delayBins(end))],[{'tiff'} {'pdf'}]);
%         load('CA1_DelayTiming','allDelayTiming')
    end
    
    if doLoc
        figure
        set(gcf,'position',[500 450 400 length(labels).*300])
        locRemapVals = repmat({[]},[1 length(labels)]);
        step = 0.01;
        cm = [[[0:step:1]'; ones(length(0:step:1),1)] ...
            [[0:step:1 1:-step:0]'] ...
            flipud([[0:step:1]'; ones(length(0:step:1),1)])];
        for i = 1:length(labels)
            subplot(length(labels),1,i)
            plotVals = cat(3,crossLocComp{i,:});
            mask = nansum(~isnan(plotVals),3)<10;
            plotVals(repmat(mask,[1 1 length(plotVals(1,1,:))])) = nan;
            imagesc(nanmedian(plotVals,3))
            alpha(double(~isnan(nanmedian(plotVals,3))))
            colormap(cm)
            caxis([-0.15 0.15])
            colorbar
            axis equal
            axis off
            
            tmp = reshape(plotVals,[numel(plotVals(:,:,1)) length(plotVals(1,1,:))]);
            locRemapVals{i} = tmp;
            apv = nan(length(tmp(:,1)),1);
            for j = 1:length(tmp(:,1))
%                 [h pval ci tstat] = ttest(tmp(j,~isnan(tmp(j,:)))',0);
                if all(isnan(tmp(j,:)))
                    continue
                end
                [pval h stat] = signrank(tmp(j,~isnan(tmp(j,:)))');
                apv(j) = pval;
            end
            apv = reshape(apv,size(plotVals(:,:,1)));

            for j = 1:length(plotVals(1,:,1))
                for k = 1:length(plotVals(:,1,1))
                    if apv(j,k) < 0.05
                        text(k,j+0.075,'*','horizontalalignment','center',...
                            'verticalalignment','middle','fontname','arial',...
                            'fontweight','bold','fontsize',9)
                    end
                end
            end
        end
        drawnow
        saveFig(gcf,[root '/RemappingXLocation'],[{'tiff'} {'pdf'}]);
    end
    
    if doLoc && doBehavioralDiagnostics
        figure
        set(gcf,'position',[500 450 200 length(labels).*200])
        for i = 1:length(labels)
            subplot(1,length(labels),i)
            tmpA = abs(nanmedian(locBiasVals{i},2));
            tmpB = nanmedian(locRemapVals{i},2);
            mask = nansum(~isnan(locRemapVals{i}),2)<10;
            tmpA(mask) = [];
            tmpB(mask) = [];
            scatter(tmpA(~isnan(tmpA)&~isnan(tmpB)), ...
                tmpB(~isnan(tmpA)&~isnan(tmpB)))
            set(gca,'xlim',[0 8])
            lsline
            [r pval] = corr(tmpA(~isnan(tmpA)&~isnan(tmpB)), ...
                tmpB(~isnan(tmpA)&~isnan(tmpB)));
            title(['r=' num2str(r) ', p=' num2str(pval)]);
            axis square
        end
        drawnow
        saveFig(gcf,[root '/RemappingXLocation_VS_Velocity'],[{'tiff'} {'pdf'}]);
%         save([root '/rrVSloc']);
%         load(['Plots/Summary/TwoSmall_CA3/rrVSloc']);
    end
    
%     figure
%     set(gcf,'position',[500 450 500 length(labels).*200])
%     for i = 1:length(labels)
%         pv = nan(25,4);
%         for j = 1:length(upiece)
%             pv(1:length(crossAllComp{i,j}(:,1)),j) = -diff(crossAllComp{i,j},[],2);
%         end
%         tmp1 = pv(1:7,:);
%         tmp2 = pv(8:14,:);
%         tmp3 = pv(15:end,:);
%         mkGraph([{tmp1(:)} {tmp2(:)} {tmp3(:)}])
%         
%         [pval h stat] = signrank(tmp1(~isnan(tmp1)))
%         [pval h stat] = signrank(tmp2(~isnan(tmp2)))
%         [pval h stat] = signrank(tmp3(~isnan(tmp3)))
%         
%         [pval h stat] = ranksum(tmp1(~isnan(tmp1)),tmp2(~isnan(tmp2)))
%         [pval h stat] = ranksum(tmp3(~isnan(tmp3)),tmp2(~isnan(tmp2)))
%         [pval h stat] = ranksum(tmp1(~isnan(tmp1)),tmp3(~isnan(tmp3)))
%     end

    
    outP = ['Stats/' paths{1}(slashInds(1):slashInds(2)-1) '.txt'];
    checkP(outP);
    fid = fopen(outP,'w');
    fprintf(fid,'\n\n\t\t\tPop Vec Remapping t-test\n');
    allGroupComps = [];
    for group = 1:length(labels)
        tmp = cat(1,crossAllComp{group,:});
        fprintf(fid,['\nWithin Group:  ' labels{group}]);
        [h pval ci tstat] = ttest(tmp(:,1),tmp(:,2));
        fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.12f ',tstat.df,tstat.tstat,pval);
        allGroupComps = [allGroupComps {diff(fliplr(tmp),[],2)}];
    end
    for i = 1:length(allGroupComps)
        for j = i+1:length(allGroupComps)
            fprintf(fid,['\nAcross Group:  ' labels{i} ' vs. ' labels{j}]);
            [h pval ci tstat] = ttest2(allGroupComps{i},allGroupComps{j});
            fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.12f ',tstat.df,tstat.tstat,pval);
        end
    end
    
    fprintf(fid,'\n\n\t\t\tPop Vec Remapping Nonparamentrics\n');
    allGroupComps = [];
    for group = 1:length(labels)
        tmp = cat(1,crossAllComp{group,:});
        fprintf(fid,['\nWithin Group (sign-rank):  ' labels{group}]);
        [pval h stat] = signrank(tmp(:,1),tmp(:,2));
        if ~isfield(stat,'zval')
                stat.zval = nan;
            end
        fprintf(fid,'\n\tSigned-rank(%i)=%0.3f, p=%0.12f ',stat.signedrank,stat.zval,pval);
        allGroupComps = [allGroupComps {diff(fliplr(tmp),[],2)}];
    end
    for i = 1:length(allGroupComps)
        for j = i+1:length(allGroupComps)
            fprintf(fid,['\nAcross Group (rank-sum):  ' labels{i} ' vs. ' labels{j}]);
            [pval h stat] = ranksum(allGroupComps{i},allGroupComps{j});
            if ~isfield(stat,'zval')
                stat.zval = nan;
            end
            fprintf(fid,'\n\tSigned-rank(%i)=%0.3f, p=%0.12f ',stat.ranksum,stat.zval,pval);
        end
    end
    
    tmp = cat(1,amfr{i,:});
    fprintf(fid,'\n\n\t\t\tCell-wise rate change t-test\n');
    allGroupComps = [];
    for group = 1:length(labels)
        tmp = cat(1,amfr{group,:});
        fprintf(fid,['\nWithin Group:  ' labels{group}]);
        [h pval ci tstat] = ttest(tmp(:,1),tmp(:,2));
        fprintf(fid,'\n\tt(%i)=%0.3f, p=%0.12f ',tstat.df,tstat.tstat,pval);
        allGroupComps = [allGroupComps {diff(fliplr(tmp),[],2)}];
    end
    
    tmp = cat(1,amfr{i,:});
    fprintf(fid,'\n\n\t\t\tCell-wise rate change Nonparametric\n');
    allGroupComps = [];
    for group = 1:length(labels)
        tmp = cat(1,amfr{group,:});
        fprintf(fid,['\nWithin Group (sign-rank):  ' labels{group}]);
        [pval h stat] = signrank(tmp(:,1),tmp(:,2));
        fprintf(fid,'\n\tSigned-rank(%i)=%0.3f, p=%0.12f ',stat.signedrank,stat.zval,pval);
        allGroupComps = [allGroupComps {diff(fliplr(tmp),[],2)}];
    end
    
%     if length(labels) > 1
%         fprintf(fid,'\n\n\t\t\tCell-wise other coding prop differences\n');
%         figure
%         set(gcf,'position',[50 50 1000 500])
%         subplot(2,4,1)
%         cumHist([{cat(1,tmfr{1,:})} {cat(1,tmfr{2,:})}]);
%         title('Mean Firing Rate')
%         [pval h stat] = ranksum(cat(1,tmfr{1,:}),cat(1,tmfr{2,:}));
%         fprintf(fid,'\n\tMFR Signed-rank(%i)=%0.3f, p=%0.12f ',stat.ranksum,stat.zval,pval);
%         subplot(2,4,5)
%         cumHist([{cat(1,tpfr{1,:})} {cat(1,tpfr{2,:})}]);
%         title('Peak Firing Rate')
%         [pval h stat] = ranksum(cat(1,tpfr{1,:}),cat(1,tpfr{2,:}));
%         fprintf(fid,'\n\tPFR Signed-rank(%i)=%0.3f, p=%0.12f ',stat.ranksum,stat.zval,pval);
%         subplot(2,4,2)
%         cumHist([{cat(1,tsi{1,:})} {cat(1,tsi{2,:})}]);
%         title('All Spatial Information')
%         [pval h stat] = ranksum(cat(1,tsi{1,:}),cat(1,tsi{2,:}));
%         fprintf(fid,'\n\tAll SI Signed-rank(%i)=%0.3f, p=%0.12f ',stat.ranksum,stat.zval,pval);
%         subplot(2,4,6)
%         cumHist([{cat(1,tsip{1,:})} {cat(1,tsip{2,:})}]);
%         title('Place Cell Spatial Information')
%         [pval h stat] = ranksum(cat(1,tsip{1,:}),cat(1,tsip{2,:}));
%         fprintf(fid,'\n\tPlace SI Signed-rank(%i)=%0.3f, p=%0.12f ',stat.ranksum,stat.zval,pval);
%         drawnow
%         subplot(2,4,3)
%         cumHist([{cat(1,tcorr{1,:})} {cat(1,tcorr{2,:})}],[-1:0.01:1]);
%         title('All Split Half')
%         [pval h stat] = ranksum(cat(1,tcorr{1,:}),cat(1,tcorr{2,:}));
%         fprintf(fid,'\n\tAll split-half Signed-rank(%i)=%0.3f, p=%0.12f ',stat.ranksum,stat.zval,pval);
%         subplot(2,4,7)
%         cumHist([{cat(1,tcorrp{1,:})} {cat(1,tcorrp{2,:})}],[-1:0.01:1]);
%         title('Place Cell Split Half')
%         [pval h stat] = ranksum(cat(1,tsip{1,:}),cat(1,tsip{2,:}));
%         fprintf(fid,'\n\tPlace split-half Signed-rank(%i)=%0.3f, p=%0.12f ',stat.ranksum,stat.zval,pval);
%         drawnow
%         subplot(2,4,4)
%         cumHist([{cat(1,tcorrpval{1,:})} {cat(1,tcorrpval{2,:})}],[0:0.01:1]);
%         title('All Split Half')
%         [pval h stat] = ranksum(cat(1,tcorrpval{1,:}),cat(1,tcorrpval{2,:}));
%         fprintf(fid,'\n\tAll split-half Signed-rank(%i)=%0.3f, p=%0.12f ',stat.ranksum,stat.zval,pval);
%         subplot(2,4,8)
%         cumHist([{cat(1,tcorrpvalp{1,:})} {cat(1,tcorrpvalp{2,:})}],[0:0.01:1]);
%         title('Place Cell Split Half')
%         [pval h stat] = ranksum(cat(1,tcorrpvalp{1,:}),tcorrpvalp(1,tsip{2,:}));
%         fprintf(fid,'\n\tPlace split-half Signed-rank(%i)=%0.3f, p=%0.12f ',stat.ranksum,stat.zval,pval);
%         drawnow
%         saveFig(gcf,[root '/Additional_Coding_Properties'],[{'tiff'} {'pdf'}]);
%     end
    
    
    fclose(fid);
    
    out.labels = labels;
    out.amfr = amfr;
    out.icorr = icorr;
    out.afw = afw;
    out.cellCounts = cellCounts;
    out.crossAllComp = crossAllComp;
    out.crossAllCompHallway = crossAllCompHallway;
    out.crossAllCompDelayed = crossAllCompDelayed;
    
    dataP = ['AggregatedData' paths{1}(slashInds(1):slashInds(2)-1)];
    checkP(dataP);
    save(dataP);
    
    
%     signal = (diff(cat(1,amfr{:}),[],2));
%     tmp = cat(1,tsnr{:});
%     figure(1)
%     set(gcf,'position',[50 50 800 400])
%     subplot(2,4,1)
%     scatter(tmp(:,1),signal)
%     corr(tmp(:,1),signal)
%     lsline
%     xlabel('SNR')
%     ylabel('MFR Difference')
%     subplot(2,4,2)
%     scatter([tmp(:,3)-tmp(:,2)]./nanmax(tmp(:,3),tmp(:,2)),signal)
%     corr([tmp(:,3)-tmp(:,2)]./nanmax(tmp(:,3),tmp(:,2)),signal)
%     lsline
%     xlabel('Change in SNR (%)')
%     ylabel('MFR Difference')
%     subplot(2,4,3)
%     mkGraph(tmp,[{'Whole Trial'} {'First Half'} {'Second Half'}])
%     ylabel('SNR')
%     subplot(2,4,4)
%     mkGraph([tmp(:,3)-tmp(:,2)]./nanmax(tmp(:,3),tmp(:,2)))
%     ylabel('Change in SNR (%)')
%     set(gca,'ylim',[-0.2 0.2])
%     
%     tmp = cat(1,tsip{:});
%     subplot(2,4,5)
%     scatter(tmp(:,1),signal)
%     corr(tmp(:,1),signal)
%     lsline
%     xlabel('SI p-value')
%     ylabel('MFR Difference')
%     
%     drawnow
%     saveFig(gcf,[root '/Revision/MFRvsProps'],[{'tiff'} {'pdf'}]);
%     
%     
%     
%     splitters = cat(1,aggVals{:});
%     signal = (cat(1,amfr{:}));
%     clear binSplit
%     binSplit(:,1) = floor(splitters(:,4).*20)+1;
%     binSplit(:,2) = floor(splitters(:,5).*20)+1;
%     binnedSNR = nan(nanmax(binSplit(:)));
%     figure
%     set(gcf,'position',[50 50 800 400])
%     subplot(2,4,1)
%     hist3(splitters(:,4:5),'nbins',[20 20]);
%     subplot(2,4,4)
%     cumHist([{signal(splitters(:,1)>10,1)} {signal(splitters(:,1)>10,2)}],[0:0.01:1])
%     subplot(2,4,5)
%     cumHist([{signal(:,1)} {signal(:,2)}],[0:0.01:1])
%     subplot(2,4,6)
%     cumHist([{signal(splitters(:,4)<0.05,1)} {signal(splitters(:,4)<0.05,2)}],[0:0.01:1])
%     subplot(2,4,7)
%     cumHist([{signal(splitters(:,5)<0.05,1)} {signal(splitters(:,5)<0.05,2)}],[0:0.01:1])
%     subplot(2,4,8)
%     cumHist([{signal(splitters(:,4)<0.05 & splitters(:,5)<0.05,1)} ...
%         {signal(splitters(:,4)<0.05 & splitters(:,5)<0.05,2)}],[0:0.01:1])
%     
end

function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end

function help_plotMaps(totalMaps,p)

    doK = [8 4];

    for part = 0:floor(length(totalMaps(1,1,:,1))/prod(doK))

        figure(1)
        set(gcf,'position',[50 50 900 1350])
        for k = 1:prod(doK)
            if part.*prod(doK)+k > length(totalMaps(1,1,:,1))
                break
            end
%             tmp = [totalMaps(:,:,part.*prod(doK)+k,1) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,2)];

            tmp = [totalMaps(:,:,part.*prod(doK)+k,1) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,3); ...
                nan(1,length(totalMaps(1,:,1,1)).*2+1); ...
                totalMaps(:,:,part.*prod(doK)+k,2) nan(size(totalMaps,1),1) totalMaps(:,:,part.*prod(doK)+k,4)];

            subplot(doK(1),doK(2),k)

            imagesc(tmp)
            colormap jet
            caxis([0 nanmax(tmp(:))])
            alpha(double(~isnan(tmp)))
            axis equal
            axis off    
        end

        slashInds = find(ismember(p{1},'/'));
        outP = ['Plots/DifferentiatedMaps/' p{1}(slashInds+1:end-4) '_Partition_' num2str(part+1)];
        saveFig(gcf,outP,'tiff')
        saveFig(gcf,outP,'pdf')
        close all
        drawnow
    end
end
























