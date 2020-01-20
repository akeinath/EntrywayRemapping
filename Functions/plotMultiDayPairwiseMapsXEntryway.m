function plotMultiDayPairwiseMapsXEntryway(paths)
    
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
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);
        s = load(paths{isM(1)});
        alignMap = s.alignment.alignmentMap;
        am = [];
        isPC = [];
        for si = 1:length(sessions)
            s = load(sessions{si});  
            [isIn isMostRecent indexSinceIn] = isInROI(s.processed.p,s.processed.roi.door);
            [isInRoom] = isInROI(s.processed.p,s.processed.roi.room);
            allMasks = [{isMostRecent(1,isInRoom)}  {isMostRecent(2,isInRoom)}];
            gT = s.processed.trace;
            mapA = mkTraceMaps(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks{1},[13 13]);
            mapB = mkTraceMaps(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks{2},[13 13]);
            map2 = cat(1,mapA,nan([1 length(mapA(1,:,1)) length(mapA(1,1,:))]),mapB);
            norm = nanmax(nanmax(map2,[],1),[],2);
            map2 = map2./repmat(norm,size(map2(:,:,1)));
            am = [am {cat(2,nan(length(map2(:,1,1)),1,length(map2(1,1,:))),map2)}];
            isPC = [isPC {s.processed.splithalf.roomXdoors.p<0.05}];
        end
        
        for si = 1:length(alignMap)-1
            for sj = si+1 %:length(alignMap)
        
                isGood = isPC{si}(alignMap{si,sj}(:,1)) | isPC{sj}(alignMap{si,sj}(:,2));
                
                m = cat(2,am{si}(:,:,alignMap{si,sj}(isGood,1)), ...
                    am{sj}(:,:,alignMap{si,sj}(isGood,2)));
                
                doK = [8 8];
                for part = 0:floor(length(m)/prod(doK))

                    figure(1)
                    set(gcf,'position',[50 50 900 900],'color','k')
                    for k = 1:prod(doK)
                        if part.*prod(doK)+k > length(m(1,1,:))
                            break
                        end
                        subplot(doK(1),doK(2),k)
                        imagesc(m(:,:,part.*prod(doK)+k))
                        colormap('jet')
                        alpha(double(~isnan(m(:,:,part.*prod(doK)+k))))
        %                 title(num2str(nanmax(nanmax(m(:,:,part.*prod(doK)+k)))))
        %                 drawnow
                        set(gca,'ydir','normal')
                        axis equal
                        axis off
                    end

                    slashInds1 = find(ismember(sessions{si},'/'));
                    slashInds2 = find(ismember(sessions{sj},'/'));
                    outP = ['Plots/PairwiseAlignedCellMapsXEntryway/' upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end) '/' ...
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