function quantBehavior(paths,doPlot)

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
    agg = [];
    
    %%%%%% ENTRIES NEED TO LAST AT LEAST 1sec to count
    aggLabels = [{'Mouse'} {'Group'} {'% Time in Compartment'} ...
        {'Number of Entries'} {'Median Duration per Entry (s)'} {'Entryway Bias'} ...
        {'Median Cumulative Distance per Entry (cm)'} ...
        {'Data included per half and entryway post-sampling matching (s)'} ...
        {'Area sampled per half and entryway post-sampling matching (%)'}];
    for mi = 1:length(upiece)
        fprintf(['\n\n\tMouse:  ' num2str(upiece{mi}) '\n']) 
        isM = find(ismember(piece,upiece(mi)));
        magg = [];
        for si = 1:length(isM);
            fprintf(['\n\t\tSession:  ' paths{isM(si)}])
            s = load(paths{isM(si)},'processed','properties');
%             s.processed.p = imfilter(s.processed.p,fspecial('gauss',[1 15],2),'same','replicate');
            %id group
            group = find(ismember(labels,paths{isM(si)}(find(ismember(paths{isM(si)},'_'),1,'last')+1:end-4)));
            
            if isempty(group)
                group = 3;
            end
            
            [isIn isMostRecent] = isInROI(s.processed.p,s.processed.roi.door);
            [isInRoom blah blah2 blah3 indexSinceIn distanceSinceIn] = ...
                isInROI(s.processed.p,s.processed.roi.room);
            indexSinceIn = indexSinceIn.*(1./30);
            half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;
            
            numberOfEntries = nansum(indexSinceIn==1);
            
            tmp = indexSinceIn(isInRoom);
            durs = tmp(logical([diff(tmp)<0 1]));
            durs(durs<1) = [];
            medianDuration = nanmedian(durs);
            
            entrywayBias = nanmax(nanmean(isMostRecent(1,indexSinceIn==1)),nanmean(isMostRecent(2,indexSinceIn==1))) ./ ...
                nanmin(nanmean(isMostRecent(1,indexSinceIn==1)),nanmean(isMostRecent(2,indexSinceIn==1)));
            
            tmp = indexSinceIn(isInRoom);
            tmp2 = distanceSinceIn(isInRoom);
            doEnds = find(logical([diff(tmp)<0 1]));
            durs = tmp(doEnds);
            doEnds(durs<1) = [];
            medianDistance = nanmedian(tmp2(doEnds));
            
            %%% Quantify coverage
            if isfield(s.processed,'exclude')
                inds = (s.processed.splithalf.roomXdoors.p<=0.05 & s.processed.exclude.SFPs);
            else
                inds = (s.processed.splithalf.roomXdoors.p<=0.05);
            end

            gT = s.processed.trace(inds,:);

            [isIn isMostRecent] = isInROI(s.processed.p,s.processed.roi.door);
            [isInRoom] = ...
                isInROI(s.processed.p,s.processed.roi.room);
            half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;

            allMasks = repmat({[]},[1 4]);
            for i = 1:2
                allMasks{i} = [isMostRecent(i,isInRoom) & half(1,isInRoom)];
                allMasks{i+2} = [isMostRecent(i,isInRoom) & ~half(1,isInRoom)];
            end

            queryMask = false(4);
            queryMask(1:2,3:4) = true;

            [a samp c d e f g m] = getMatchedMapsNMasks(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks,queryMask);

            sampOverlap = nansum(samp(:)>0)./12.^2;
            
            slashInds = find(ismember(upiece{mi},'/'));
            magg = [magg; {upiece{mi}(slashInds(end)+1:end)} {s.properties.trial(1:8)} labels(group) ...
                {nanmean(isInRoom).*100} {numberOfEntries} {medianDuration} {entrywayBias} ...
                {medianDistance} {nansum(samp(:))./30} {sampOverlap.*100}];
        end
        dg = unique(magg(:,3));
        for gi = 1:length(dg)
            dat = reshape(cat(1,magg{ismember(magg(:,3),dg(gi)),4:end}), ...
                [nansum(ismember(magg(:,3),dg(gi))) length(magg(1,:))-3]);
            tmp = [nanmean(dat,1); nanstd(dat,1)];
            
            t2 = [];
            for j = 1:length(tmp(1,:))
                t2 = [t2 {sprintf([char(181) '=%0.2f ' char(963) '=%0.2f'],[tmp(1,j) tmp(2,j)])}];
            end
            
            agg = [agg; magg(1,1)  dg(gi) t2];
        end
    end
    xlswrite('BehavioralQuantification',[aggLabels; agg])
end


























