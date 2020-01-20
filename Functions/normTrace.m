function normTrace(paths)
    clc
    fprintf(['\t\nNorming traces with binarizing.\n\n'])
    warning off all
    for p = paths'
        s = load(p{1});
        fprintf(['\t\t' num2str(p{1}) '\n\t\t\tModeling spike trains:      '])
        

        t = s.calcium.FiltTraces(s.processed.validTraceFrames(:,1),:)';
        isGood = find(diff(t(1,:),[],2)~=0);
        for k = 1:length(t(:,1))
            t(k,isGood(1):isGood(end)) = linterp(isGood,t(k,isGood),isGood(1):isGood(end));
        end
        
%         zt = bsxfun(@rdivide,bsxfun(@minus,t,nanmean(t,2)),nanstd(t,[],2));
%         dzt = [nan(length(zt(:,1)),1) diff(zt,[],2)];
%         isRise = dzt > 0 & zt > 1;
%         s.processed.trace = isRise;
%         s.processed.traceModel = {'Binarized'};

        t = detrend(t')'; % norm to the median
        
        dzt = [nan(length(t(:,1)),1) diff(t,[],2)];
        
        isRise = false(size(t));
        for k = 1:length(t(:,1))
            bsd = nanstd(t(k,t(k,:)<0))./sqrt(1-2./pi);
            zi = t(k,:)./bsd;
            isRise(k,:) = zi>2 & dzt(k,:)>0;
        end
        s.processed.trace = isRise;
        s.processed.traceModel = {'HalfNorm_Binarized'};
        
        
%         t = t';
%         spikes = nan(size(s.calcium.FiltTraces(s.processed.validTraceFrames(:,1),:)));
%         for i = 1:length(s.calcium.FiltTraces(1,:))
%             str_f = sprintf('%6.1f',100.*i/length(s.calcium.FiltTraces(1,:)));
%             fprintf([repmat('\b',[1 6]) str_f])
%             [blah spikes(:,i)] = deconvolveCa(t(:,i),'ar2');
%         end
%         fprintf('\n')
%         s.processed.trace = spikes';
%         s.processed.traceModel = {'AutoReg2'};
        
        save(p{1},'-struct','s','-v7.3');
    end
end