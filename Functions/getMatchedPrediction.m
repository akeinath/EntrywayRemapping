function [m os val ival locval locival] = getMatchedMapsNMasks(p,ts,masks)
    
    for k = 1:length(masks)
        [m(:,:,:,k) s(:,:,k)] = mkTraceMaps(p,ts,masks{k},[]); % [15 15]
    end
    
    os = min(s,[],3);
    
    nsim = 100;
    if all(os(:)==0)
        val = nan;
    	return
    end
    
    for k = 1:length(masks)
        om(:,:,:,:,k) = mkTraceMaps(p,ts,masks{k},[],os,nsim); % [15 15]
    end

    stepSize = 30;
    prediction = nan(1,floor(length(p)./stepSize));
    for sti = 1:stepSize:length(p)
        
        for k = 1:length(masks)
            [m(:,:,:,k) s(:,:,k)] = mkTraceMaps(p,ts,masks{k},[]); % [15 15]
        end
        os = min(s,[],3);
        for k = 1:length(masks)
            om(:,:,:,:,k) = mkTraceMaps(p(:,sMask),ts(:,sMask),masks{k}(:,sMask),[],os,nsim); % [15 15]
        end
    end
    
end

























