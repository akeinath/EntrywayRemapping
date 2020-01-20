function plotMultiDayPairwiseMaps(paths)
    
    warning off all
    
    clc
    fprintf('\n')
    %%% Reliability constaint
    
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
    
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end)) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);
        s = load(paths{isM(1)});
        alignMap = s.alignment.alignmentMap;
        am = repmat({[]},[1 length(sessions)]);
        isPC = repmat({[]},[1 length(sessions)]);
        for si = 1:length(sessions)
            s = load(sessions{si},'processed');
            slashInds = find(ismember(sessions{si},'/'));
            fprintf(['\t\tPreloading:  ' sessions{si}(slashInds(end)+1:end-4) '\n'])
            
            gT = s.processed.trace;
            am{si} = mkTraceMaps(s.processed.p,gT);
            isPC{si} = s.processed.splithalf.roomXdoors.p <= 0.05;
        end  
        
        for si = 1:length(sessions)
            for sj = si+1:length(sessions)
                tmp1 = am{si};
                tmp2 = am{sj};
                
                if isempty(alignMap{si,sj})
                    continue
                end
                
                
                isGood = isPC{si}(alignMap{si,sj}(:,1)) | isPC{sj}(alignMap{si,sj}(:,2));
                toPlot = [tmp1(:,:,alignMap{si,sj}(isGood,1)) nan(length(tmp1(:,1,1)),4,nansum(isGood)) ...
                    tmp2(:,:,alignMap{si,sj}(isGood,2))];
                
                doK = [8 4];

                for part = 0:floor(length(toPlot(1,1,:,1))/prod(doK))

                    figure(1)
                    set(gcf,'position',[50 50 900 1350])
                    for k = 1:prod(doK)
                        if part.*prod(doK)+k > length(toPlot(1,1,:))
                            break
                        end
                        
                        subplot(doK(1),doK(2),k)
                        imagesc(toPlot(:,:,part.*prod(doK)+k))
                        colormap jet
%                         caxis([0 nanmax(nanmax(toPlot(:,:,part.*prod(doK)+k)))])
                        alpha(double(~isnan(toPlot(:,:,part.*prod(doK)+k))))
                        axis equal
                        axis off    
                    end
                    slashInds1 = find(ismember(sessions{si},'/'));
                    slashInds2 = find(ismember(sessions{sj},'/'));
                    outP = ['Plots/PairwiseAlignedCellMaps/' upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end) '/' ...
                        [sessions{si}(slashInds1(end)+1:end-4) '_vs_' sessions{sj}(slashInds2(end)+1:end-4)] ...
                        '_Part_' num2str(part)];
                    saveFig(gcf,outP,[{'tiff'} {'pdf'}])
                    close all
                    drawnow
                end
            end
        end
    end
end