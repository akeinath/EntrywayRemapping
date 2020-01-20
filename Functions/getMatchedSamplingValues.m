function [val ival mfr pfr pfshift fieldCenters bmrs] = getMatchedSamplingValues(p,ts,masks,queryMask)

    if nargin < 4 || isempty(queryMask)
        queryMask = true(length(masks));
    end
    ival = nan;
    mfr = nan;
    pfr = [];
    pfshift = [];
    bmrs = [];
    for k = 1:length(masks)
        [m(:,:,:,k) s(:,:,k)] = mkTraceMaps(p,ts,masks{k});
    end
    
    os = min(s,[],3);
    
    
    nsim = 100;
    if all(os(:)==0)
        val = nan;
    	return
    end
    
    for k = 1:length(masks)
        [om(:,:,:,:,k) uns(:,:,:,:,k) unm(:,:,:,:,k)] = mkTraceMaps(p,ts,masks{k},[15 25],os,nsim); % [14 14]
    end
    
%     unm = om;

    if length(masks) == 4
        unm(:,:,:,:,1:2) = bsxfun(@rdivide,unm(:,:,:,:,1:2),nanmax(nansum(nansum(unm(:,:,:,:,1:2),1),2),[],5)); % normalize
        unm(:,:,:,:,3:4) = bsxfun(@rdivide,unm(:,:,:,:,3:4),nanmax(nansum(nansum(unm(:,:,:,:,3:4),1),2),[],5)); % normalize
    end
    if length(masks) == 2
        unm(:,:,:,:,1:2) = bsxfun(@rdivide,unm(:,:,:,:,1:2),nanmax(nansum(nansum(unm(:,:,:,:,1:2),1),2),[],5)); % normalize
    end
    
%     unm = bsxfun(@rdivide,unm,nanmax(nanmax(nanmax(unm,[],1),[],2),[],5)); % normalize
    
%     ros = repmat(os,[1 1 length(om(1,1,:,1,1)) nsim]);

    ival = nan([length(masks) length(masks) length(ts(:,1))]);
    mfr = nan([length(masks) length(masks) length(ts(:,1))]);
    pfr = nan([length(masks) length(masks) length(ts(:,1))]);
    bmrs = nan([length(masks) length(masks) 4]);
    pfshift = nan([length(masks) length(masks) length(ts(:,1))]);
    fieldCenters = nan([length(ts(:,1)) 2 length(masks)]);
    for ki = 1:length(masks)
        ma = om(:,:,:,:,ki);
        
        normA = nansum(nansum(ma,1),2);
        normA(normA==0) = 1;
        tma = bsxfun(@rdivide,ma,normA);
        [x y] = meshgrid(1:length(ma(1,:,1,1)),1:length(ma(:,1,1,1)));
        ax = permute(nansum(nansum(tma.*repmat(x,[1 1 length(tma(1,1,:,1)) ...
            length(tma(1,1,1,:))]),1),2),[3 4 1 2]);
        ay = permute(nansum(nansum(tma.*repmat(y,[1 1 length(tma(1,1,:,1)) ...
            length(tma(1,1,1,:))]),1),2),[3 4 1 2]);
        fieldCenters(:,1,ki) = nanmean(ax,2);
        fieldCenters(:,2,ki) = nanmean(ay,2);
        for kj = ki+1:length(masks)
            if ~queryMask(ki,kj)
                continue
            end
            mb = om(:,:,:,:,kj);
            
            %%% norm to max to eliminate relative weight
%             normA = nanmax(nanmax(ma,[],1),[],2);
%             normA(normA==0) = 1;
%             normB = nanmax(nanmax(mb,[],1),[],2);
%             normB(normB==0) = 1;
% 
%             ma = bsxfun(@rdivide,ma,normA);
%             mb = bsxfun(@rdivide,mb,normB);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            a = reshape(ma,[numel(ma(:,:,:,1)) nsim]);
            b = reshape(mb,[numel(ma(:,:,:,1)) nsim]);
            a(isnan(a(:,1)),:) = [];
            b(isnan(b(:,1)),:) = [];
            tmp = corr(a,b);
            val(ki,kj) = nanmean(tmp(logical(eye(length(tmp)))));
            
            if nargout > 1
                for k = 1:length(ma(1,1,:,1))
                    a = reshape(ma(:,:,k,:),[numel(ma(:,:,1,1)) nsim]);
                    b = reshape(mb(:,:,k,:),[numel(ma(:,:,1,1)) nsim]);
                    a(isnan(a(:,1)),:) = [];
                    b(isnan(b(:,1)),:) = [];
                    tmp = corr(a,b);
                    ival(ki,kj,k) = nanmean(tmp(logical(eye(length(tmp)))));
                end
            end
            
            if nargout > 2
                sa = nansum(nansum(unm(:,:,:,:,ki),1),2);
                sb = nansum(nansum(unm(:,:,:,:,kj),1),2);
                if length(masks) == 2
                    mfr(ki,kj,:) = nanmean(sa-sb,4);
                else
                    mfr(ki,kj,:) = nanmean(abs(sa-sb),4);
                end
            end
            
            if nargout > 3
                sa = nanmax(nanmax(unm(:,:,:,:,ki),[],1),[],2);
                sb = nanmax(nanmax(unm(:,:,:,:,kj),[],1),[],2);
                pfr(ki,kj,:) = nanmean(abs(sa-sb),4);
            end
            
            if nargout > 4
                
                normA = nansum(nansum(ma,1),2);
                normA(normA==0) = 1;
                normB = nansum(nansum(mb,1),2);
                normB(normB==0) = 1;

                tma = bsxfun(@rdivide,ma,normA);
                tmb = bsxfun(@rdivide,mb,normB);
                
                [x y] = meshgrid(1:length(ma(1,:,1,1)),1:length(ma(:,1,1,1)));
                
                ax = permute(nansum(nansum(tma.*repmat(x,[1 1 length(tma(1,1,:,1)) ...
                    length(tma(1,1,1,:))]),1),2),[3 4 1 2]);
                ay = permute(nansum(nansum(tma.*repmat(y,[1 1 length(tma(1,1,:,1)) ...
                    length(tma(1,1,1,:))]),1),2),[3 4 1 2]);
                bx = permute(nansum(nansum(tmb.*repmat(x,[1 1 length(tma(1,1,:,1)) ...
                    length(tma(1,1,1,:))]),1),2),[3 4 1 2]);
                by = permute(nansum(nansum(tmb.*repmat(y,[1 1 length(tma(1,1,:,1)) ...
                    length(tma(1,1,1,:))]),1),2),[3 4 1 2]);
                
                pfshift(ki,kj,:) = nanmean(sqrt((ax-bx).^2 + (ay-by).^2),2);
                
                
%                 a = reshape(ma,[numel(ma(:,:,1,1)) length(ma(1,1,:,1)) nsim]);
%                 b = reshape(mb,[numel(ma(:,:,1,1)) length(ma(1,1,:,1)) nsim]);
%                 
%                 [blah apeak] = nanmax(a,[],1);
%                 apeak = permute(apeak,[3 2 1]);
%                 [xa ya] = ind2sub(size(ma(:,:,1,1)),apeak);
%                 
%                 [blah bpeak] = nanmax(b,[],1);
%                 bpeak = permute(bpeak,[3 2 1]);
%                 [xb yb] = ind2sub(size(ma(:,:,1,1)),bpeak);
%                 
%                 pfshift(ki,kj,:) = nanmean(sqrt((xa-xb).^2 + (ya-yb).^2),1);
            end
            
            if nargout > 6
                rval = nan(1,4);
                for rot = 0:3
                    rmb = imrotate(mb,rot.*90,'crop');
                    a = reshape(ma,[numel(ma(:,:,:,1)) nsim]);
                    b = reshape(rmb,[numel(ma(:,:,:,1)) nsim]);
                    a(isnan(a(:,1)),:) = [];
                    b(isnan(b(:,1)),:) = [];
                    tmp = corr(a,b);
                    rval(rot+1) = nanmean(tmp(logical(eye(length(tmp)))));
                end
                [blah ind] = nanmax(rval);
%                 bmrs(ki,kj) = [ind-1].*90;
                bmrs(ki,kj,:) = rval;
            end
        end
    end
end
