function [unmatch match] = getMatchedMapsBehavioralBiases(p,masks,queryMask)
    if nargin < 3 || isempty(queryMask)
        queryMask = true(length(masks));
    end

    for k = 1:length(masks)
        [unmatch(:,:,:,k) s(:,:,k)] = mkBehavMaps(p,masks{k},[15 15]);
    end
    
    os = min(s,[],3);
    
    nsim = 100;
    if all(os(:)==0)
        val = nan;
    	return
    end
    
    for k = 1:length(masks)
        match(:,:,:,:,k) = mkBehavMaps(p,masks{k},[15 15],os,nsim);
    end
    match = permute(nanmean(match,4),[1 2 3 5 4]);
end

























