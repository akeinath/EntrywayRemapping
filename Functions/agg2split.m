function agg2split(folder)
    if nargin<1 || isempty(folder)
        folder = 'AggregatedData';
    end
    p = getFilePaths(folder,'.mat');

    s{1} = load(p{2});
    s{2} = load(p{4});
    s{3} = load(p{1});
    s{4} = load(p{3});
    
    for si = 1:length(s)
        
        
        close all
        drawnow
        
        figure
        set(gcf,'position',[50 50 1800 300.*length(s{si}.labels)])
        for gi = 1:length(s{si}.labels)
            signal = (cat(1,s{si}.amfr{gi,:}));
            splitters = cat(1,s{si}.aggVals{gi,:});

            clear binSplit
            binSplit(:,1) = floor(splitters(:,4).*20)+1;
            binSplit(:,2) = floor(splitters(:,5).*20)+1;
            
            subplot(length(s{si}.labels),6,[(gi-1).*6]+[1:2])
            N = hist3(splitters(:,4:5),'nbins',[20 20]);
            imagesc(N./nansum(N(:)))
            colormap('hot')
            colorbar
            set(gca,'ydir','normal','xtick',[],'ytick',[])
            xlabel('Split-half correlation (p-val)')
            ylabel('Spatial information content (SIC; p-val)')
            axis equal
            axis off
            
            subplot(length(s{si}.labels),6,(gi-1).*6+3)
            h = cumHist([{signal(:,1)} {signal(:,2)}],[0:0.01:1]);
            xlabel(sprintf('Split-half change\nin mean firing rate (%%)'),'interpreter','none')
            ylabel('Cumulative Proportion')
            set(gca,'xticklabel',100.*cellfun(@str2num,get(gca,'xticklabel')))
            title('All Cells')
            [pval h stat] = signrank(signal(:,1),signal(:,2));
%             text(0.95,0.12,sprintf('Z=%0.2f\np=%0.12f',abs(stat.zval),pval),'fontname','arial',...
%                 'fontsize',9,'fontweight','normal','horizontalalignment','right');
            
            subplot(length(s{si}.labels),6,(gi-1).*6+4)
            cumHist([{signal(splitters(:,4)<0.05,1)} {signal(splitters(:,4)<0.05,2)}],[0:0.01:1]);
            xlabel(sprintf('Split-half change\nin mean firing rate (%%)'),'interpreter','none')
            ylabel('Cumulative Proportion')
            set(gca,'xticklabel',100.*cellfun(@str2num,get(gca,'xticklabel')))
            title('Split-half corr p < 0.05')
            [pval h stat] = signrank(signal(splitters(:,4)<0.05,1),signal(splitters(:,4)<0.05,2));
%             text(0.95,0.12,sprintf('Z=%0.2f\np=%0.12f',abs(stat.zval),pval),'fontname','arial',...
%                 'fontsize',9,'fontweight','normal','horizontalalignment','right');
            
            
            subplot(length(s{si}.labels),6,(gi-1).*6+5)
            cumHist([{signal(splitters(:,5)<0.05,1)} {signal(splitters(:,5)<0.05,2)}],[0:0.01:1]);
            xlabel(sprintf('Split-half change\nin mean firing rate (%%)'),'interpreter','none')
            ylabel('Cumulative Proportion')
            set(gca,'xticklabel',100.*cellfun(@str2num,get(gca,'xticklabel')))
            title('SIC p < 0.05')
            [pval h stat] = signrank(signal(splitters(:,5)<0.05,1),signal(splitters(:,5)<0.05,2));
%             text(0.95,0.12,sprintf('Z=%0.2f\np=%0.12f',abs(stat.zval),pval),'fontname','arial',...
%                 'fontsize',9,'fontweight','normal','horizontalalignment','right');
            
            subplot(length(s{si}.labels),6,(gi-1).*6+6)
            h = cumHist([{signal(splitters(:,4)<0.05 & splitters(:,5)<0.05,1)} ...
                {signal(splitters(:,4)<0.05 & splitters(:,5)<0.05,2)}],[0:0.01:1]);          
            xlabel(sprintf('Split-half change\nin mean firing rate (%%)'),'interpreter','none')
            ylabel('Cumulative Proportion')
            set(gca,'xticklabel',100.*cellfun(@str2num,get(gca,'xticklabel')))
            title('Split-half corr p < 0.05 & SIC p < 0.05')
            legend([h{1}(1) h{2}(1)],[{'Same'} ...
                {'Diff.'} {'Shuffled'}],'location','east','fontname','arial',...
                'fontsize',9,'fontweight','bold','color','none','box','off')
            [pval h stat] = signrank(signal(splitters(:,4)<0.05 & splitters(:,5)<0.05,1), ...
                signal(splitters(:,4)<0.05 & splitters(:,5)<0.05,2));
%             text(0.95,0.12,sprintf('Z=%0.2f\np=%0.12f',abs(stat.zval),pval),'fontname','arial',...
%                 'fontsize',9,'fontweight','normal','horizontalalignment','right');
        end
        
        drawnow
        saveFig(gcf,['Plots/Revision/SplitAggregated_Exp_' num2str(si)],[{'tiff'} {'pdf'}]);
        
    end
end