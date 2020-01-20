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
    pThresh = 0.05; % 0.05;
    crossAllComp = repmat({[]},[length(labels) length(upiece) 4]);
    for mi = 1:length(upiece)
        fprintf(['\n\n\tMouse:  ' num2str(upiece{mi}) '\n']) 
        isM = find(ismember(piece,upiece(mi)));
        for si = 1:length(isM);
            fprintf(['\n\t\tSession:  ' paths{isM(si)}])
            s = load(paths{isM(si)},'processed','properties');
            
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
            doInclude = true(length(s.processed.trace(:,1)),1);
            if isfield(s.processed,'exclude')
                gT = s.processed.trace(s.processed.exclude.SFPs,:);
                doInclude = s.processed.exclude.SFPs;
            else
                gT = s.processed.trace;
            end
                                  
            % s.processed.splithalf.roomXdoors_si.p<=pThresh & s.processed.splithalf.roomXdoors.p<=pThresh
            
            [isIn isMostRecent] = isInROI(s.processed.p,s.processed.roi.door);
            [isInRoom blah blah2 blah3 indexSinceIn distanceSinceIn] = ...
                isInROI(s.processed.p,s.processed.roi.room);
            indexSinceIn = indexSinceIn.*(1./30);
            half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;

            allMasks = repmat({[]},[1 4]);
            for i = 1:2
                allMasks{i} = [isMostRecent(i,isInRoom) & half(1,isInRoom)];
                allMasks{i+2} = [isMostRecent(i,isInRoom) & ~half(1,isInRoom)];
            end
            
            queryMask = false(4);
            queryMask(1:2,3:4) = true;
          
            [map samp allComp] = getMatchedMapsNMasks(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks,queryMask);

            crossAllComp{group,mi,1} = cat(1,crossAllComp{group,mi,1},...
                [help_getMaskedVals(allComp(1:2,3:4,:),[true false; false true]) ...
                help_getMaskedVals(allComp(1:2,3:4,:),[false true; true false])]);
            
            [map samp allComp] = getMatchedMapsNMasks(s.processed.p(:,isInRoom), ...
                gT(s.processed.splithalf.roomXdoors.p(doInclude)<=pThresh,isInRoom),allMasks,queryMask);

            crossAllComp{group,mi,2} = cat(1,crossAllComp{group,mi,2},...
                [help_getMaskedVals(allComp(1:2,3:4,:),[true false; false true]) ...
                help_getMaskedVals(allComp(1:2,3:4,:),[false true; true false])]);
            
            [map samp allComp] = getMatchedMapsNMasks(s.processed.p(:,isInRoom), ...
                gT(s.processed.splithalf.wholemap_si.p(doInclude)<=pThresh,isInRoom),allMasks,queryMask);

            crossAllComp{group,mi,3} = cat(1,crossAllComp{group,mi,3},...
                [help_getMaskedVals(allComp(1:2,3:4,:),[true false; false true]) ...
                help_getMaskedVals(allComp(1:2,3:4,:),[false true; true false])]);
            
            [map samp allComp] = getMatchedMapsNMasks(s.processed.p(:,isInRoom), ...
                gT(s.processed.splithalf.wholemap_si.p(doInclude)<=pThresh & ...
                s.processed.splithalf.roomXdoors.p(doInclude)<=pThresh,isInRoom),allMasks,queryMask);

            crossAllComp{group,mi,4} = cat(1,crossAllComp{group,mi,4},...
                [help_getMaskedVals(allComp(1:2,3:4,:),[true false; false true]) ...
                help_getMaskedVals(allComp(1:2,3:4,:),[false true; true false])]);

        end
    end
    close all    
    
    slashInds = find(ismember(paths{1},'/'));
    root = ['Plots/REVISION/Summary' paths{1}(slashInds(1):slashInds(2)-1)];
    
    eliminate = all(cellfun(@isempty,crossAllComp),2);
    crossAllComp(eliminate(:,:,1),:,:) = [];
    labels(eliminate(:,:,1)) = [];
    
    toPlotDiff = repmat({[]},[1 length(crossAllComp(1,1,:)).*length(crossAllComp(:,1,1))]);
    toPlot = repmat({[]},[2 length(crossAllComp(1,1,:)).*length(crossAllComp(:,1,1))]);
    for mi = 1:length(crossAllComp(1,1,:))
        for gi = 1:length(crossAllComp(:,1,1))
            tmp = cat(1,crossAllComp{gi,:,mi});
            toPlot{1,(gi-1).*length(crossAllComp(1,1,:))+mi} = tmp(:,1);
            toPlot{2,(gi-1).*length(crossAllComp(1,1,:))+mi} = tmp(:,2);
            toPlotDiff{(gi-1).*length(crossAllComp(1,1,:))+mi} = tmp(:,1)-tmp(:,2);
        end
    end
    
    slashInds = find(ismember(paths{1},'/'));
    root = ['Plots/REVISION/MultiSelection_PV/' paths{1}(slashInds(1):slashInds(2)-1)];
    
    figure
    set(gcf,'position',[50 50 350.*length(labels) 250])
    mkGraph(toPlot,[],[{'Same'} {'Diff'}])
    set(gca,'ylim',[0 1])
    ylabel('Correlation (r)')
    saveFig(gcf,[root '_Total'],[{'pdf'} {'tiff'}]); 
    
    figure
    set(gcf,'position',[50 50 175.*length(labels) 250])
    mkGraph(toPlotDiff,[])
    ylabel('[PV_s - PV_d]')
    saveFig(gcf,[root '_Diff'],[{'pdf'} {'tiff'}]); 
    
    tmp = toPlotDiff;
    apv = [];
    for i = 1:length(tmp(1,:))
        [pval h stat] = signrank(tmp{i});
        apv = [apv [pval; 0]];
    end
    
    outP = ['Stats/MultiCriteria' paths{1}(slashInds(1):slashInds(2)-1) '.txt'];
    checkP(outP);
    fid = fopen(outP,'w');
    allGroupComps = [];
    fprintf(fid,'\n\n\t\t\tPop Vec Remapping Nonparamentrics Multi Crit\n');
    fprintf(fid,'\n\tSigned-rank Zs > %0.3f, p < %0.12f ',nanmin(apv(2,:)),nanmax(apv(1,:)));
    
    if length(tmp(1,:))>4
        apv = nan([length(tmp(1,:)) length(tmp(1,:)) 2]);
        for i = 1:length(tmp(1,:))
            for j = i+1:length(tmp(1,:))
                [pval h stat] = ranksum(tmp{i},tmp{j});
                apv(i,j,:) = [pval 0];
            end
        end
        tmp2 = apv(1:4,5:8,:);
        
    end
    fclose(fid);
    
    a = load('a');
    a = a.tmp;
    b = tmp;
    apv = nan(2,8);
    for i = 1:length(tmp(1,:))
        [pval h stat] = ranksum(a{i},b{i});
        apv(:,i) = [pval stat.zval];
    end
%     save('a')
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
























