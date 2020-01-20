function selectionCriteria2(paths,varargin)
    clc
    close all
    drawnow
    
    if isempty(varargin)
        varargin = [{'Wholemap'}];
    end
    varargin = cellfun(@lower,varargin,'uniformoutput',false);
    
    warning off all
    if isempty(gcp)
        parpool('local',7);
    end
    pctRunOnAll warning off all
    
    
    fprintf(['Computing split-half reliability:\n']);
    tmp = [repmat({'\n\t'},[1 length(varargin)]); varargin];
    fprintf(cat(2,tmp{:},'\n\n'))
    velThresh = -2;
    nsims = 500; %500
    minShift = 900;
    for p = paths'
        s = load(p{1});
        didChange = false;
%         s.processed.p = imfilter(s.processed.p,fspecial('gauss',[1 15],2),'same','replicate');
%         s.processed = rmfield(s.processed,'splithalf');
        fprintf(['\t' num2str(p{1}) '\n'])    
        if isfield(s.processed,'exclude')
            include = s.processed.exclude.SFPs;
        else
            include = true(length(s.processed.trace(:,1)),1);
        end
        if ismember({'wholemap'},varargin)
            %%% In-room split half reliability, ignore doorway

            
%             if isfield(s.processed,'splithalf') && isfield(s.processed.splithalf,'wholemap')
%                 continue
%             end
            didChange = true;
            
            half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;
            allMasks = [{half} {~half}];
            [a b c ival] = getMatchedMapsNMasks(s.processed.p,s.processed.trace,allMasks);
            s.processed.splithalf.wholemap.val = ival;

            null = nan(length(s.processed.trace(:,1)),nsims);
            fprintf('\t\tWhole Map, computing null...  ')
            tic
            P = s.processed.p;
            T = s.processed.trace;
            parfor sim = 1:nsims
                gT = circshift(T,[0 minShift+randi(length(T(1,:))-minShift.*2)]);

                [map samp allComp ival] = getMatchedMapsNMasks(P,gT,allMasks);
                null(:,sim) = permute(ival(1,2,:),[3 2 1]);
            end
            durat = toc;
            fprintf([num2str(durat) ' sec']);
            s.processed.splithalf.wholemap.p = 1-nanmean(bsxfun(@gt,...
                permute(s.processed.splithalf.wholemap.val(1,2,:),[3 2 1]),[null(repmat(include,[1 nsims]))]'),2);


            figure(1)
            set(gcf,'position',[50 50 250 250])
            cumHist(s.processed.splithalf.wholemap.p,[0:0.01:1]);
            hold on
            plot([0 1],[0 1],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',1)
            ylabel('Cumulative Proportion')
            xlabel('P-Value')
            drawnow
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/SplitHalf/WholeMap/' p{1}(slashInds+1:end-4)];
            saveFig(gcf,outP,'tiff')
            close all
            drawnow

            fprintf('\n\t\tReliable spatial cells:  %i of %i (%0.2f%%)\n',[nansum(s.processed.splithalf.wholemap.p <= 0.05) ...
                length(s.processed.splithalf.wholemap.val(1,2,:)) nanmean(s.processed.splithalf.wholemap.p <= 0.05).*100]);
        end

        if ismember({'room'},varargin)
            %%% In-room split half reliability, ignore doorway

            
            didChange = true;
            isInRoom = isInROI(s.processed.p,s.processed.roi.room);

            half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;
            allMasks = [{half(:,isInRoom)} {~half(:,isInRoom)}];
            [a b c ival] = getMatchedMapsNMasks(s.processed.p(:,isInRoom),s.processed.trace(:,isInRoom),allMasks);
            s.processed.splithalf.room.val = ival;


            null = nan(length(s.processed.trace(:,1)),nsims);
            fprintf('\t\tWithin room, ignoring doorways, computing null...  ')
            tic
            P = s.processed.p(:,isInRoom);
            T = s.processed.trace(:,isInRoom);
            parfor sim = 1:nsims
                gT = circshift(T,[0 minShift+randi(length(T(1,:))-minShift.*2)]);

                [map samp allComp ival] = getMatchedMapsNMasks(P,gT,allMasks);

                null(:,sim) = permute(ival(1,2,:),[3 2 1]);
            end
            durat = toc;
            fprintf([num2str(durat) ' sec']);
            s.processed.splithalf.room.p = 1-nanmean(bsxfun(@gt,...
                permute(s.processed.splithalf.room.val(1,2,:),[3 2 1]),[null(repmat(include,[1 nsims]))]'),2);


            figure(1)
            set(gcf,'position',[50 50 250 250])
            cumHist(s.processed.splithalf.room.p,[0:0.01:1]);
            hold on
            plot([0 1],[0 1],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',1)
            ylabel('Cumulative Proportion')
            xlabel('P-Value')
            drawnow
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/SplitHalf/Room/' p{1}(slashInds+1:end-4)];
            saveFig(gcf,outP,'tiff')
            close all
            drawnow

            fprintf('\n\t\tReliable spatial cells:  %i of %i (%0.2f%%)\n',[nansum(s.processed.splithalf.room.p <= 0.05) ...
                length(s.processed.splithalf.room.val(1,2,:)) nanmean(s.processed.splithalf.room.p <= 0.05).*100]);
        end
        
        if ismember({'roomxdoor'},varargin)
            %%% In-room split half reliability, include doorway

            vel = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
            vel = imfilter(vel,fspecial('gauss',[1 30],10),'same','replicate');
            
            didChange = true;
            [isIn isMostRecent indexSinceIn] = isInROI(s.processed.p,s.processed.roi.door);
            isInRoom = isInROI(s.processed.p,s.processed.roi.room);

            half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;
            allMasks = repmat({[]},[1 length(s.processed.roi.door(1,:)).*2]);
            for doorA = 1:length(s.processed.roi.door(1,:))
                allMasks{doorA} = [isMostRecent(doorA,isInRoom) & half(1,isInRoom) & vel(1,isInRoom)>velThresh];
                allMasks{doorA+length(s.processed.roi.door(1,:))} = [isMostRecent(doorA,isInRoom) & ~half(1,isInRoom) & vel(1,isInRoom)>velThresh];
                
%                 h2 = cumsum(isMostRecent(doorA,isInRoom));
%                 allMasks{doorA} = [isMostRecent(doorA,isInRoom) & h2<h2(end)./2];
%                 allMasks{doorA+length(s.processed.roi.door(1,:))} = [isMostRecent(doorA,isInRoom) & ~(h2<h2(end)./2)];

            end
            [a b c ivals] = getMatchedMapsNMasks(s.processed.p(:,isInRoom),s.processed.trace(:,isInRoom),allMasks);

            tmp = ivals(1:end/2,end/2+1:end,:);
            tmp = permute(tmp,[3 2 1]);
            allComp = reshape(tmp,length(tmp(:,1,1)),[]);

            allComp = [nanmean(allComp(:,[1 4]),2) nanmean(allComp(:,[2 3]),2)];

            allComp = sort(allComp,2,'descend');
            actual = allComp;
            s.processed.splithalf.roomXdoors.vals = allComp;


            null = nan(length(s.processed.trace(:,1)),length(s.processed.roi.door(1,:)),nsims);
            fprintf('\t\tWithin room by doorways, computing null...  ')
            tic
            P = s.processed.p(:,isInRoom);
            T = s.processed.trace(:,isInRoom);
            clear ivals
            parfor sim = 1:nsims
                gT = circshift(T,[0 minShift+randi(length(T(1,:))-minShift.*2)]);

                [map samp allComp ivals] = getMatchedMapsNMasks(P,gT,allMasks);

                tmp = ivals(1:end/2,end/2+1:end,:);
                tmp = permute(tmp,[3 2 1]);
                allComp = reshape(tmp,length(tmp(:,1,1)),[]);

                allComp = [nanmean(allComp(:,[1 4]),2) nanmean(allComp(:,[2 3]),2)];
                null(:,:,sim) = allComp;
            end
            durat = toc;
            fprintf([num2str(durat) ' sec']);
            null = sort(null,2,'descend');
            s.processed.splithalf.roomXdoors.null = null;    
            tmp = nanmax(null,[],2);
%             s.processed.splithalf.roomXdoors.p = 1-nanmean(bsxfun(@gt,...
%                 nanmax(s.processed.splithalf.roomXdoors.vals,[],2)',tmp(:)))';
            s.processed.splithalf.roomXdoors.p = 1-nanmean(bsxfun(@gt,...
                nanmax(s.processed.splithalf.roomXdoors.vals,[],2)',permute(tmp,[3 1 2])))';

            figure(1)
            set(gcf,'position',[50 50 250 250])
            cumHist(s.processed.splithalf.roomXdoors.p,[0:0.01:1]);
            hold on
            plot([0 1],[0 1],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',1)
            ylabel('Cumulative Proportion')
            xlabel('P-Value')
            drawnow
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/SplitHalf/RoomXDoorway/' p{1}(slashInds+1:end-4)];
            saveFig(gcf,outP,'tiff')
            close all
            drawnow

            fprintf('\n\t\tReliable spatial cells:  %i of %i (%0.2f%%)\n',[nansum(s.processed.splithalf.roomXdoors.p <= 0.05) ...
                length(s.processed.splithalf.roomXdoors.vals(:,1)) nanmean(s.processed.splithalf.roomXdoors.p <= 0.05).*100]);
        end
        
        if ismember({'hallwayxdoor'},varargin)
            %%% In-hallway split half reliability, include doorway

            didChange = true;
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
%                 if oneEntrance(start)~=oneEntrance(start+stop)
                    if oneEntrance(start)
                        isDir(1,start:start+stop) = true;
                    else
                        isDir(2,start:start+stop) = true;
                    end
%                 end
                isGood(start:start+stop) = false;
            end

            half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;

            allMasks = [{isDir(1,isIn) & half(isIn)} {isDir(2,isIn)& half(isIn)} ...
                {isDir(1,isIn) & ~half(isIn)} {isDir(2,isIn)& ~half(isIn)}];
            [maps samp allComp ivals] = getMatchedLinMapsNMasks(cp(isIn),s.processed.trace(:,isIn),allMasks);
            
            tmp = ivals(1:end/2,end/2+1:end,:);
            tmp = permute(tmp,[3 2 1]);
            allComp = reshape(tmp,length(tmp(:,1,1)),[]);


            allComp = sort(allComp,2,'descend');
            actual = allComp;
            s.processed.splithalf.hallwayXdoors.vals = allComp;
            
            null = nan(length(s.processed.trace(:,1)),length(s.processed.roi.door(1,:)).^2,nsims);
            fprintf('\t\tWithin hallway by doorways, computing null...  ')
            tic
            P = cp(:,isIn);
            T = s.processed.trace(:,isIn);
            clear ivals
            parfor sim = 1:nsims
                gT = circshift(T,[0 minShift+randi(length(T(1,:))-minShift.*2)]);

                [maps samp allComp ivals] = getMatchedLinMapsNMasks(P,gT,allMasks);

                tmp = ivals(1:end/2,end/2+1:end,:);
                tmp = permute(tmp,[3 2 1]);
                allComp = reshape(tmp,length(tmp(:,1,1)),[]);

                null(:,:,sim) = allComp;
            end
            durat = toc;
            fprintf([num2str(durat) ' sec']);
            null = sort(null,2,'descend');
            s.processed.splithalf.hallwayXdoors.null = null;    
            tmp = nanmax(null,[],2);
%             s.processed.splithalf.hallwayXdoors.p = 1-nanmean(bsxfun(@gt,...
%                 nanmax(s.processed.splithalf.hallwayXdoors.vals,[],2)',tmp(:)))';
            s.processed.splithalf.hallwayXdoors.p = 1-nanmean(bsxfun(@gt,...
                nanmax(s.processed.splithalf.hallwayXdoors.vals,[],2)',permute(tmp,[3 1 2])))';


            figure(1)
            set(gcf,'position',[50 50 250 250])
            cumHist(s.processed.splithalf.hallwayXdoors.p,[0:0.01:1]);
            hold on
            plot([0 1],[0 1],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',1)
            ylabel('Cumulative Proportion')
            xlabel('P-Value')
            drawnow
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/SplitHalf/HallwayXDoorway/' p{1}(slashInds+1:end-4)];
            saveFig(gcf,outP,'tiff')
            close all
            drawnow

            fprintf('\n\t\tReliable spatial cells:  %i of %i (%0.2f%%)\n',[nansum(s.processed.splithalf.hallwayXdoors.p <= 0.05) ...
                length(s.processed.splithalf.hallwayXdoors.vals(:,1)) nanmean(s.processed.splithalf.hallwayXdoors.p <= 0.05).*100]);
        end
        
        if ismember({'rate_roomxdoor'},varargin)
            
            didChange = true;
            %%% In-room split half reliability, include doorway

            [isIn isMostRecent indexSinceIn] = isInROI(s.processed.p,s.processed.roi.door);
            isInRoom = isInROI(s.processed.p,s.processed.roi.room);
            half = 1:length(s.processed.p(1,:)) < length(s.processed.p(1,:))./2;
            gT = s.processed.trace(:,isInRoom);

            allMasks = repmat({[]},[1 4]);
            for i = 1:2
                allMasks{i} = [isMostRecent(i,isInRoom) & half(1,isInRoom)];
                allMasks{i+2} = [isMostRecent(i,isInRoom) & ~half(1,isInRoom)];
            end

            
            [blah1 ivals mfr pfr] = getMatchedSamplingValues(s.processed.p(:,isInRoom),gT,...
                [{isMostRecent(1,isInRoom)} {isMostRecent(2,isInRoom)}]);
            s.processed.splithalf.rate.room.vals = permute(mfr(1,2,:),[3 1 2]);
            
            null = nan(length(s.processed.trace(:,1)),nsims);
            fprintf('\t\tRate within room by doorways, computing null...  ')
            tic
            T = s.processed.trace(:,isInRoom);
            clear ivals
            for sim = 1:nsims
                gT = circshift(T,[0 minShift+randi(length(T(1,:))-minShift.*2)]);

                
                [blah1 ivals mfr pfr] = getMatchedSamplingValues(s.processed.p(:,isInRoom),gT,...
                    [{isMostRecent(1,isInRoom)} {isMostRecent(2,isInRoom)}]);
                null(:,sim) = permute(mfr(1,2,:),[3 1 2]);
            end
            durat = toc;
            fprintf([num2str(durat) ' sec']);
            s.processed.splithalf.rate.room.null = null;    
            tmp = nanmax(null,[],2);
            s.processed.splithalf.rate.room.p = 1-nanmean(bsxfun(@gt,...
                s.processed.splithalf.rate.room.vals,null),2);

            figure(1)
            set(gcf,'position',[50 50 250 250])
            cumHist(s.processed.splithalf.rate.room.p,[0:0.01:1]);
            hold on
            plot([0 1],[0 1],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',1)
            ylabel('Cumulative Proportion')
            xlabel('P-Value')
            drawnow
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/SplitHalf/Rate_InRoom/' p{1}(slashInds+1:end-4)];
            saveFig(gcf,outP,'tiff')
            close all
            drawnow

            fprintf('\n\t\tReliable spatial cells:  %i of %i (%0.2f%%)\n\t\t\tPercent of place cells:  %0.2f%%\n',[nansum(s.processed.splithalf.rate.room.p <= 0.05) ...
                length(s.processed.splithalf.rate.room.vals(:,1)) nanmean(s.processed.splithalf.rate.room.p <= 0.05).*100 ...
                nanmean(s.processed.splithalf.rate.room.p(s.processed.splithalf.roomXdoors.p<=0.05)<=0.05).*100]);
        end
        
        if ismember({'wholemap_si'},varargin)
            %%% In-room split half reliability, include doorway

%             if ~(isfield(s.processed,'splithalf') && isfield(s.processed.splithalf,'wholemap_si'))
                didChange = true;
            
                isInRoom = isInROI(s.processed.p,s.processed.roi.room);

                actual = [];
                [blah b c] = mkTraceMaps(s.processed.p(:,isInRoom),...
                    s.processed.trace(:,isInRoom));
                b = b./nansum(b(:));
                actual = permute(nansum(nansum(bsxfun(@times,b,(c./repmat(nanmean(nanmean(c,1),2),[size(b)])) .* ...
                    log(c./repmat(nanmean(nanmean(c,1),2),[size(b)]))),1),2),[3 1 2]);

                fprintf('\t\tWhole map (spatial information), computing null...  ')
                tic
                P = s.processed.p(:,isInRoom);
                T = s.processed.trace(:,isInRoom);
                clear ivals
                null = nan([length(actual) nsims]);
                parfor sim = 1:nsims
                    gT = circshift(T,[0 minShift+randi(length(T(1,:))-minShift.*2)]);

                    [blah b c] = mkTraceMaps(P,gT);
                    b = b./nansum(b(:));
                    si = permute(nansum(nansum(bsxfun(@times,b,(c./repmat(nanmean(nanmean(c,1),2),[size(b)])) .* ...
                        log(c./repmat(nanmean(nanmean(c,1),2),[size(b)]))),1),2),[3 1 2]);
                    null(:,sim) = si;
                end
                durat = toc;
                s.processed.splithalf.wholemap_si.vals = actual;
                s.processed.splithalf.wholemap_si.p = ...
                    1 - nanmean(bsxfun(@gt,actual,null),2);

                figure(1)
                set(gcf,'position',[50 50 250 250])
                cumHist(s.processed.splithalf.wholemap_si.p,[0:0.01:1]);
                hold on
                plot([0 1],[0 1],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',1)
                ylabel('Cumulative Proportion')
                xlabel('P-Value')
                drawnow
                slashInds = find(ismember(p{1},'/'));
                outP = ['Plots/SplitHalf/WholeMap_SI/' p{1}(slashInds+1:end-4)];
                saveFig(gcf,outP,'tiff')
                close all
                drawnow


                fprintf('\n\t\tReliable spatial cells:  %i of %i (%0.2f%%)\n',[nansum(s.processed.splithalf.wholemap_si.p <= 0.05) ...
                    length(s.processed.splithalf.wholemap_si.p(:,1)) nanmean(s.processed.splithalf.wholemap_si.p <= 0.05).*100]);            
%             end
        end
        
        if ismember({'roomxdoor_si'},varargin)
            %%% In-room split half reliability, include doorway

            didChange = true;
            [isIn isMostRecent indexSinceIn] = isInROI(s.processed.p,s.processed.roi.door);
            isInRoom = isInROI(s.processed.p,s.processed.roi.room);
            
            actual = [];
            for k = 1:2
                [blah b c] = mkTraceMaps(s.processed.p(:,isInRoom),...
                    s.processed.trace(:,isInRoom),isMostRecent(k,isInRoom));
                b = b./nansum(b(:));
                si = permute(nansum(nansum(bsxfun(@times,b,(c./repmat(nanmean(nanmean(c,1),2),[size(b)])) .* ...
                    log(c./repmat(nanmean(nanmean(c,1),2),[size(b)]))),1),2),[3 1 2]);
                actual(:,k) = si;
            end
            
            fprintf('\t\tWithin hallway by doorways (spatial information), computing null...  ')
            tic
            P = s.processed.p(:,isInRoom);
            T = s.processed.trace(:,isInRoom);
            clear ivals
            null = nan([size(actual) nsims]);
            mask = isMostRecent(:,isInRoom);
            parfor sim = 1:nsims
                gT = circshift(T,[0 minShift+randi(length(T(1,:))-minShift.*2)]);

                for k = 1:2
                    [blah b c] = mkTraceMaps(P,...
                        gT,mask(k,:));
                    b = b./nansum(b(:));
                    si = permute(nansum(nansum(bsxfun(@times,b,(c./repmat(nanmean(nanmean(c,1),2),[size(b)])) .* ...
                        log(c./repmat(nanmean(nanmean(c,1),2),[size(b)]))),1),2),[3 1 2]);
                    null(:,k,sim) = si;
                end
            end
            durat = toc;
            s.processed.splithalf.roomXdoors_si.vals = actual;
            s.processed.splithalf.roomXdoors_si.p = ...
                1 - nanmean(bsxfun(@gt,nanmax(actual,[],2),permute(nanmax(null,[],2),[1 3 2])),2);
            
            figure(1)
            set(gcf,'position',[50 50 250 250])
            cumHist(s.processed.splithalf.roomXdoors_si.p,[0:0.01:1]);
            hold on
            plot([0 1],[0 1],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',1)
            ylabel('Cumulative Proportion')
            xlabel('P-Value')
            drawnow
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/SplitHalf/RoomXDoors_SI/' p{1}(slashInds+1:end-4)];
            saveFig(gcf,outP,'tiff')
            close all
            drawnow

            
            fprintf('\n\t\tReliable spatial cells:  %i of %i (%0.2f%%)\n',[nansum(s.processed.splithalf.roomXdoors_si.p <= 0.05) ...
                length(s.processed.splithalf.roomXdoors_si.p(:,1)) nanmean(s.processed.splithalf.roomXdoors_si.p <= 0.05).*100]);            
        end

        if didChange
            save(p{1},'-struct','s','-v7.3');
        end
    end
end

function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end

function sh = help_splithalf(p,t)

    v = [0 sqrt(sum(diff(p,[],2).^2))].*30;
    m1 = mkTraceMaps(p,t,1:length(p(1,:))<length(p(1,:))./2);
    m2 = mkTraceMaps(p,t,1:length(p(1,:))>=length(p(1,:))./2);

    rm1 = nan(length(m1(1,:,1)).*length(m1(:,1,1)),length(m1(1,1,:)));
    rm2 = nan(length(m1(1,:,1)).*length(m1(:,1,1)),length(m1(1,1,:)));
    for k = 1:length(m1(1,1,:))
        tmp = m1(:,:,k);
        rm1(:,k) = tmp(:);
        tmp = m2(:,:,k);
        rm2(:,k) = tmp(:);
    end

    isBad = isnan(rm1) | isnan(rm2);

    m1t = reshape(rm1(~isBad),[nansum(~isBad(:,1)) length(m1(1,1,:))]);
    m2t = reshape(rm2(~isBad),[nansum(~isBad(:,1)) length(m1(1,1,:))]);

    xc = (corr(m1t,m2t));
    sh = xc(logical(eye(size(xc))));
end