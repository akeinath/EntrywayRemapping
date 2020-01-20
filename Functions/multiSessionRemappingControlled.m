
function multiSessionRemappingControlled(paths)
    
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
        
        allRR = repmat({[]},[length(isM) length(isM)]);
        allPC = repmat({[]},[1 length(length(isM))]);
        groups = nan(1,length(isM));
        for di = 1:1:length(isM);
            s1 = load(paths{isM(di)});
            if di==1
                alignMap = s1.alignment.alignmentMap;
            end
            inds1 = (s1.processed.splithalf.roomXdoors.p<=pThresh);
            if isfield(s1.processed,'exclude')
                inds1 = inds1 & s1.processed.exclude.SFPs;
            end
            group = find(ismember(labels,paths{isM(di)}( ...
                find(ismember(paths{isM(di)},'_'),1,'last')+1:end-4)));
            
            if isempty(group)
                group = 3;
            end
            groups(di) = group;
            
            for dj = di+1:1:length(isM);
                if isempty(alignMap{di,dj})
                    continue
                end
                s2 = load(paths{isM(dj)});
                
                if (length(alignMap{di,dj})./nanmax(length(s1.processed.trace(:,1)),...
                        length(s2.processed.trace(:,1)))) < 0.5 
                    continue
                end
                fprintf(['\t\tSession:  ' num2str(di) '\t' num2str(dj) '\n'])
                
                isInRoom1 = isInROI(s1.processed.p,s1.processed.roi.room);
                isInRoom2 = isInROI(s2.processed.p,s2.processed.roi.room);
                
                [isIn mrd1] = isInROI(s1.processed.p,(s1.processed.roi.door));
                [isIn mrd2] = isInROI(s2.processed.p,(s2.processed.roi.door));
                
                inds2 = (s2.processed.splithalf.roomXdoors.p<=pThresh);
                if isfield(s2.processed,'exclude')
                    inds2 = inds2 & s2.processed.exclude.SFPs;
                end
                
                aIsPC = inds1(alignMap{di,dj}(:,1));
                bIsPC = inds2(alignMap{di,dj}(:,2));
                
                P1 = s1.processed.p(:,isInRoom1);
                P2 = s2.processed.p(:,isInRoom2);
                T1 = s1.processed.trace(alignMap{di,dj}(:,1),isInRoom1);
                T2 = s2.processed.trace(alignMap{di,dj}(:,2),isInRoom2);
                P = [P1 P2];
                T = [T1(aIsPC|bIsPC,:) T2(aIsPC|bIsPC,:)];
                
                allMasks = repmat({[]},[1 length(s1.processed.roi.door(1,:)).*2]);
                for doorA = 1:length(s1.processed.roi.door(1,:))
                    allMasks{doorA} = [mrd1(doorA,isInRoom1) false(1,length(mrd2(doorA,isInRoom2)))];
                    allMasks{doorA+length(s1.processed.roi.door(1,:))} = ...
                        [false(1,length(mrd1(doorA,isInRoom1))) mrd2(doorA,isInRoom2)];
                end

                [blah1 ivals mfr] = getMatchedSamplingValues(P,T,allMasks);
                allRR{di,dj} = [help_getMaskedVals(mfr(1:2,3:4,:),[true false; false true]) ...
                    help_getMaskedVals(mfr(1:2,3:4,:),[false true; true false])];
            end
        end
        a = cat(1,allRR{logical(bsxfun(@times,groups==1,groups'==1))});
        b = cat(1,allRR{logical(bsxfun(@times,groups==2,groups'==2))});
        c = cat(1,allRR{logical(bsxfun(@times,groups==3,groups'==3))});
        d = cat(1,allRR{logical(bsxfun(@times,groups==2,groups'==1))});
        e = cat(1,allRR{logical(bsxfun(@times,groups==1,groups'==2))});
        allCrossRR = [{[allCrossRR{1}; a]} ...
            {[allCrossRR{2}; b]} ...
            {[allCrossRR{3}; c]} ...
            {[allCrossRR{4}; d]} ...
            {[allCrossRR{5}; e]}];
        
%         allCrossRR = [allCrossRR; cat(1,allRR{:})];

        h = cumHist(cat(1,allRR{:}),[0:0.025:1]);
%         h = cumHist(cat(1,amfr{i,:}),[0:0.025:1]);
        xlabel('Mean Rate Change (%)')
        ylabel('Count')
        legend([h{1}(1) h{2}(1)],[{'Within'} ...
            {'Across'} {'Shuffled'}],'location','southeast','fontname','arial',...
            'fontsize',9,'fontweight','bold','color','none','box','off')
        set(gca,'xticklabel',100.*cellfun(@str2num,get(gca,'xticklabel')))
        axis square
        slashInds = find(ismember(paths{isM(di)},'/'));
        outP = ['Plots/AlignedAnalyses/' paths{isM(di)}(slashInds(2)+1:end-4)];
        saveFig(gcf,outP,[{'tiff'} {'pdf'}])
        close all
        drawnow
    end
end


function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end















