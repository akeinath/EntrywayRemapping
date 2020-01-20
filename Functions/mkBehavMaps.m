function [m s] = mkBehavMaps(p,mask,dims,downsample,nsim)
    if nargin<3 || isempty(mask)
        mask = true(length(p(1,:)));
    end
    
    if nargin<5 || isempty(downsample)
        downsample = [];
    end
    if nargin<4 || isempty(downsample) || nargin<5 || isempty(nsim)
        nsim = 1;
    end
    
    
    binsize = 2.5; % 2.5
%     kern = 4; %    4
    
    vel = [0 sqrt(sum(diff(p,[],2).^2))].*30;
    vel = imfilter(vel,fspecial('gauss',[1 30],10),'same','replicate');

    tp = diff(p,[],2);
    [theta dist] = cart2pol(tp(1,:)',tp(2,:)');
    theta = [theta(1) theta'];

    p = bsxfun(@minus,p,nanmin(p')');
    p = floor(p./binsize)+1;
    if nargin<3 || isempty(dims)
        m = nan([nanmax(p') 2 nsim]);
    else
        m = nan([dims 2 nsim]);
    end
    s = zeros(size(m(:,:,1)));
    
    p(:,~mask) = [];
    vel(:,~mask) = [];
    theta(:,~mask) = [];
    
    up = unique(p','rows');
    up(any(isnan(up),2),:) = [];
    for loc = up'
        isGood = find([loc(1)==p(1,:) & loc(2)==p(2,:)]);
        if ~isempty(downsample)
            repGood = nan(nsim,downsample(loc(1),loc(2)));
            for si = 1:nsim
                repGood(si,:) = isGood(randperm(length(isGood),downsample(loc(1),loc(2))));
            end
            repGood = repGood';
            
            tmp = reshape(vel(:,repGood),[1 downsample(loc(1),loc(2)) nsim]);
            m(loc(1),loc(2),1,:) = permute(nanmean(tmp,2),[4 2 1 3]);
            
            tmp = reshape(theta(:,repGood),[1 downsample(loc(1),loc(2)) nsim]);
            if isempty(tmp)
                continue
            end
%             m(loc(1),loc(2),2,:) = permute(circ_r(tmp,[],[],2),[4 2 1 3]);
        else
            m(loc(1),loc(2),1) = nanmean(vel(:,isGood),2);
            
%             m(loc(1),loc(2),2) = circ_r(theta(:,isGood)');
        end
        s(loc(1),loc(2)) = length(isGood);
    end
%     unm = m;
%     bad = isnan(m);
%     m(bad) = 0;    
%     m = imfilter(m,fspecial('gauss',round([kern.*5 kern.*5]./binsize),kern./binsize),'same');
%     m(bad) = nan;

    
%     m = bsxfun(@times,m,s);
%     m = imfilter(m,fspecial('gauss',round([kern.*5 kern.*5]./binsize),kern./binsize),'same');
%     s = imfilter(s,fspecial('gauss',round([kern.*5 kern.*5]./binsize),kern./binsize),'same');
%     m = bsxfun(@rdivide,m,s);
end






















