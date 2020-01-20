function plotMatchedMaps(paths)
    for p = paths'
        s = load(p{1});
        
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

        [a b c d e f g m] = getMatchedMapsNMasks(s.processed.p(:,isInRoom),gT(:,isInRoom),allMasks,queryMask);

        help_plotMaps2(m(:,:,:,:,1),p)
    end
end

function help_plotMaps2(totalMaps,p)

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
        outP = ['Plots/DifferentiatedMaps_MatchedSampling/' p{1}(slashInds+1:end-4) '_Partition_' num2str(part+1)];
        saveFig(gcf,outP,'tiff')
        saveFig(gcf,outP,'pdf')
        close all
        drawnow
    end
end