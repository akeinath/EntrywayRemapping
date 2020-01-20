
function multiSessionRemapping(paths)
    
    clc
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
    pThresh = 0.05;
    allCrossRR = [{[]} {[]} {[]} {[]} {[]}];
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}) '\n'])
        
        isM = find(ismember(piece,upiece(mi)));
        
        allRR = repmat({[]},[1 length(length(isM))]);
        allPC = repmat({[]},[1 length(length(isM))]);
        groups = nan(1,length(isM));
        for di = 1:1:length(isM);
            s = load(paths{isM(di)});
            group = find(ismember(labels,paths{isM(di)}( ...
                find(ismember(paths{isM(di)},'_'),1,'last')+1:end-4)));
            
            if isempty(group)
                group = 3;
            end
            
            groups(di) = group;
            
            if di==1
                alignMap = s.alignment.alignmentMap;
            end
            fprintf(['\t\tSession:  ' num2str(spiece{isM(di)}) '\n'])
            
            gT = s.processed.trace(:,:);
            
            [isIn isMostRecent] = isInROI(s.processed.p,s.processed.roi.door);
            [isInRoom blah blah2 blah3 indexSinceIn] = isInROI(s.processed.p,s.processed.roi.room);
            indexSinceIn = indexSinceIn.*(1./30);
            half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;

            allMasks = repmat({[]},[1 4]);
            for i = 1:2
                allMasks{i} = [isMostRecent(i,isInRoom) & half(1,isInRoom)];
                allMasks{i+2} = [isMostRecent(i,isInRoom) & ~half(1,isInRoom)];
            end
%             [blah1 ivals mfr pfr pfshift fieldCenters] = getMatchedSamplingValues(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks);
%             tmp = [help_getMaskedVals(mfr(1:2,3:4,:),[true false; false true]) ...
%                 help_getMaskedVals(mfr(1:2,3:4,:),[false true; true false])];
%             allRR{di} = tmp;

            [blah1 ivals mfr pfr pfshift fieldCenters] = getMatchedSamplingValues(...
                s.processed.p(:,isInRoom),gT(:,isInRoom),[{isMostRecent(1,isInRoom)} {isMostRecent(2,isInRoom)}]);
            allRR{di} = permute(mfr(1,2,:),[3 1 2]);
            
            inds = (s.processed.splithalf.roomXdoors.p<=pThresh);
            if isfield(s.processed,'exclude')
                inds = inds & s.processed.exclude.SFPs;
            end
            allPC{di} = inds;
        end
        
        crossRR = nan(length(isM));
        cellwiseRR = [];
        cellwiseRR_Null = [];
        nsims = 2;
        nullRR = nan([length(isM) nsims]);
        for di = 1:1:length(isM);
            for dj = di+1:1:length(isM);
                if isempty(alignMap{di,dj}) || (length(alignMap{di,dj})./...
                        nanmax(length(allRR{di}),length(allRR{dj}))) < 0.5 
                    continue
                end
% %                 a = -diff(allRR{di}(alignMap{di,dj}(:,1),:),[],2);
% %                 aIsPC = allPC{di}(alignMap{di,dj}(:,1));
% %                 b = -diff(allRR{dj}(alignMap{di,dj}(:,2),:),[],2);
% %                 bIsPC = allPC{dj}(alignMap{di,dj}(:,2));
% %                 if nansum(aIsPC|IsPC) < 30
% %                     continue
% %                 end
% %                 crossRR(di,dj) = corr(a,b,'type','Kendall');

                aIsPC = allPC{di}(alignMap{di,dj}(:,1));
                bIsPC = allPC{dj}(alignMap{di,dj}(:,2));
                a = allRR{di}(alignMap{di,dj}(:,1));
                b = allRR{dj}(alignMap{di,dj}(:,2));
%                 if nansum(aIsPC|bIsPC) < 30
%                     continue
%                 end                
                cellwiseRR = [cellwiseRR; abs(a-b)];
                crossRR(di,dj) = corr(a(aIsPC|bIsPC),b(aIsPC|bIsPC));
                
                if corr(a(aIsPC|bIsPC),b(aIsPC|bIsPC)) > 0.25
                    figure
                    set(gcf,'position',[50 50 300 300])
                    scatter(a(aIsPC|bIsPC),b(aIsPC|bIsPC),3,'k')
                    lsline
                    ylabel('MFRa - MFRb')
                    xlabel('MFRa - MFRb')
                    drawnow
                    slashInds = find(ismember(paths{isM(di)},'/'));
                    outP = ['Plots/AlignedAnalyses/RR_Corr_Examples/' paths{isM(di)}(slashInds(2)+1:end-4)];
                    slashInds = find(ismember(paths{isM(dj)},'/'));
                    outP = [outP '_' paths{isM(dj)}(slashInds(end)+1:end-4)];
                    saveFig(gcf,outP,[{'tiff'} {'pdf'}])
                    close all
                    drawnow
                end
                
                for si = 1:nsims
                    cellwiseRR_Null = [cellwiseRR_Null; abs(a-b(randperm(length(b))))];
                end
            end
        end
        
        a = crossRR(logical(bsxfun(@times,groups==1,groups'==1)));
        b = crossRR(logical(bsxfun(@times,groups==2,groups'==2)));
        c = crossRR(logical(bsxfun(@times,groups==3,groups'==3)));
        d = crossRR(logical(bsxfun(@times,groups==2,groups'==1)));
        e = crossRR(logical(bsxfun(@times,groups==1,groups'==2)));
        allCrossRR = [{[allCrossRR{1}; a(~isnan(a))]} ...
            {[allCrossRR{2}; b(~isnan(b))]} ...
            {[allCrossRR{3}; c(~isnan(c))]} ...
            {[allCrossRR{4}; d(~isnan(d))]} ...
            {[allCrossRR{5}; e(~isnan(e))]}];
        
%         close all        
%         allCrossRR = [allCrossRR; crossRR(~isnan(crossRR))];
%         
%         figure
%         set(gcf,'position',[50 50 300 300])
%         hist(crossRR(~isnan(crossRR)))
        cumHist([{cellwiseRR} {cellwiseRR_Null}])
%         drawnow
%         slashInds = find(ismember(paths{isM(1)},'/'));
%         outP = ['Plots/AlignedAnalyses/RateRemapping_Cellwise/' paths{isM(1)}(slashInds+1:end-4)];
%         saveFig(gcf,outP,[{'tiff'} {'pdf'}])
%         close all
%         drawnow
        
%         nanmean(crossRR(logical(bsxfun(@times,groups==1,groups'==1))))
%         nanmean(crossRR(logical(bsxfun(@times,groups==2,groups'==2))))
%         nanmean(crossRR(logical(bsxfun(@times,groups==3,groups'==3))))
%         
%         hist(crossRR(logical(bsxfun(@times,groups==1,groups'==1))))
%         hist(crossRR(logical(bsxfun(@times,groups==2,groups'==2))))
%         hist(crossRR(logical(bsxfun(@times,groups==3,groups'==3))))
%         
%         save(paths{isM(1)},'-struct','s','-v7.3');
    end
    figure
    set(gcf,'position',[50 50 900 200])
    subplot(1,3,1)
    hist(allCrossRR{1})
    set(gca,'xlim',[-0.5 0.5])
    subplot(1,3,2)
    hist(allCrossRR{2})
    set(gca,'xlim',[-0.5 0.5])
    subplot(1,3,3)
    hist(cat(1,allCrossRR{4:5}))
    set(gca,'xlim',[-0.5 0.5])
    drawnow;
    [h pval stat] = signrank(allCrossRR{1})
    [h pval stat] = signrank(allCrossRR{2})
    [h pval stat] = signrank(cat(1,allCrossRR{4:5}))


    figure
    set(gcf,'position',[50 50 300 200])
    hist(cat(1,allCrossRR{:}))
    set(gca,'xlim',[-0.5 0.5])
    xlabel('Correlation (r)')
    ylabel('Count')
    [h pval stat] = signrank(cat(1,allCrossRR{:}))
    
%     [h pval stat] = signrank(allCrossRR)
%     
%     figure
%     set(gcf,'position',[50 50 300 300])
%     hist(allCrossRR)
%     drawnow
end


function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end















