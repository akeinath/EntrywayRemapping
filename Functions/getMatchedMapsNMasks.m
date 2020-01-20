function [m os val ival locval locival s om] = getMatchedMapsNMasks(p,ts,masks,queryMask)
    if nargin < 4 || isempty(queryMask)
        queryMask = true(length(masks));
    end

    ival = nan(length(masks));
    locval = nan(length(masks));
    locival = nan(length(masks));
    
    for k = 1:length(masks)
        [m(:,:,:,k) s(:,:,k)] = mkTraceMaps(p,ts,masks{k},[15 15]); % [15 15]
    end
    
    os = min(s,[],3);
    
%     ival = nan([length(masks) length(masks) length(m(1,1,:,1))]);
%     val = nan(length(masks));
%     for ki = 1:length(masks)
%         for kj = ki+1:length(masks)
%             ma = m(:,:,:,ki);
%             mb = m(:,:,:,kj);
%             val(ki,kj) = corr(ma(~isnan(ma)&~isnan(mb)),...
%                 mb(~isnan(ma)&~isnan(mb)));
%             if nargout > 3
%                 
%                 mak = reshape(ma,numel(ma(:,:,1)),[]);
%                 mbk = reshape(mb,numel(mb(:,:,1)),[]);
%                 tmp = corr(mak(any(~isnan(mak),2)&any(~isnan(mbk),2),:),...
%                     mbk(any(~isnan(mak),2)&any(~isnan(mbk),2),:));
%                 ival(ki,kj,:) = (tmp(logical(eye(length(tmp)))));
%             end
%         end
%     end 

%     figure(1)
%     set(gcf,'position',[50 50 1000 200])
%     for i = 1:4
%         subplot(1,5,i)
%         imagesc(s(:,:,i))
%         colormap parula
%         alpha(double(0~=(s(:,:,i))))
%         axis equal
%         axis off    
%         caxis([0 600])
%     end
%     subplot(1,5,5)
%     imagesc(os)
%     colormap parula
%     alpha(double(0~=(os)))
%     caxis([0 600])
%     axis equal
%     axis off
%     saveFig(gcf,['Plots/SamplingExample'],[{'pdf'} {'tiff'}])
%     
%     figure(2)
%     set(gcf,'position',[50 50 400 400])
%     imagesc(os)
%     colormap parula
%     alpha(double(0~=(os)))
%     caxis([0 600])
%     axis equal
%     axis off
%     colorbar
%     saveFig(gcf,['Plots/SamplingExample_colorbar'],[{'pdf'} {'tiff'}])

    
    nsim = 100;
    if all(os(:)==0)
        val = nan;
    	return
    end
    
    for k = 1:length(masks)
        om(:,:,:,:,k) = mkTraceMaps(p,ts,masks{k},[15 15],os,nsim); % [15 15]
    end
    
%     for i = 1:3
%         normA = nanmax(nanmax(om(:,:,:,:,[i i+1]),[],1),[],2);
%         normA(normA==0) = 1;
%         om(:,:,:,:,[i i+1]) = om(:,:,:,:,[i i+1])./repmat(normA,[size(om(:,:,1,1,1)) 1 1 1]);
%     end
    
    ival = nan([length(masks) length(masks) length(ts(:,1))]);
    locval = nan([size(om(:,:,1,1,1)) length(masks) length(masks)]);
    locival = nan([size(om(:,:,:,1,1)) length(masks) length(masks)]);
    for ki = 1:length(masks)
        for kj = ki+1:length(masks)
            if ~queryMask(ki,kj)
                continue
            end

            ma = om(:,:,:,:,ki);
            mb = om(:,:,:,:,kj);
            
%             %%% norm to max within map to eliminate relative weight
%             normA = nanmax(nanmax(ma,[],1),[],2);
%             normA(normA==0) = 1;
%             normB = nanmax(nanmax(mb,[],1),[],2);
%             normB(normB==0) = 1;
% 
%             ma = bsxfun(@rdivide,ma,normA);
%             mb = bsxfun(@rdivide,mb,normB);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%             %%% norm to max across all conditions 
%             %%% to eliminate relative weight
%             
%             normer = nanmax(nanmax(nanmax(om,[],1),[],2),[],5);
%             ma = bsxfun(@rdivide,ma,normer);
%             mb = bsxfun(@rdivide,mb,normer);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            a = reshape(ma,[numel(ma(:,:,:,1)) nsim]);
            b = reshape(mb,[numel(ma(:,:,:,1)) nsim]);
            a(isnan(a(:,1)),:) = [];
            b(isnan(b(:,1)),:) = [];
            
            
            
            tmp = corr(a,b);
            val(ki,kj) = nanmedian(tmp(logical(eye(length(tmp)))));
            
%             tmp = nan(1,nsim);
%             for i = 1:nsim
%                 x = a(:,i);
%                 y = b(:,i);
%                 xy   = dot(x,y);
%                 nx   = norm(x);
%                 ny   = norm(y);
% %                 nx = 1; %%%% Just compute dot product
% %                 ny = 1;
%                 nxny = nx*ny;
%                 tmp(i)   = xy/nxny;
%             end
%             
%             val(ki,kj) = nanmean(tmp);
            
            if nargout > 3
                for k = 1:length(ma(1,1,:,1))
                    a = reshape(ma(:,:,k,:),[numel(ma(:,:,1,1)) nsim]);
                    b = reshape(mb(:,:,k,:),[numel(ma(:,:,1,1)) nsim]);
                    a(isnan(a(:,1)),:) = [];
                    b(isnan(b(:,1)),:) = [];
                    tmp = corr(a,b);
                    ival(ki,kj,k) = nanmean(tmp(logical(eye(length(tmp)))));
                end
            end
            
            if nargout > 4
                
                locmaps = nan([size(ma(:,:,1,1)) nsim]);
                for si = 1:length(ma(1,1,1,:))
                    a = reshape(ma(:,:,:,si),[numel(ma(:,:,1,1)) size(ma,3)]);
                    b = reshape(mb(:,:,:,si),[numel(ma(:,:,1,1)) size(ma,3)]);
%                     a(isnan(a(:,1)),:) = [];
%                     b(isnan(b(:,1)),:) = [];
                    tmp = corr(a',b');
                    tmp = tmp(logical(eye(length(tmp))));
                    locmaps(:,:,si) = reshape(tmp,[size(ma(:,:,1,1))]);
                end
                locval(:,:,ki,kj) = nanmean(locmaps,3);
            end
            
            if nargout > 5
                locival(:,:,:,ki,kj) = nanmean(ma - mb,4);
            end
        end
    end
end

























