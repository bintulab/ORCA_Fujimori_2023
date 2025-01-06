%%
% import /Users/tfuji/Documents/tempanalysis/EM-seq/EMseqViz_combine_all.mat
%% median of DNAmet fraction
MedDNAmet = [];
SizedFig(15,100);
visidx = [26,28,22,20  27,29,23,21];
for ii = 1:length(visidx)
    subplot(length(visidx),1,ii);
    tmp = dnamet{visidx(ii)};
    tmp(tmp == 2) = NaN;
    tmp(tmp == -1) = NaN;
    tmp = 1-nanmean(tmp,2);
    histogram(tmp,20,'BinLimits',[0,1],'normalization','probability');
    set(gca,'XTick',xax,'XTickLabel',xaxlab);
    set(gca,'XColor','k','YColor','k');
%     title(num2str(round(nanmedian(tmp),3)));
    title(dname{visidx(ii)});
    MedDNAmet = [MedDNAmet,nanmedian(tmp)];
%     title(num2str(round(sum(tmp > 0.5)/length(tmp),3)));
%     MedDNAmet = [MedDNAmet,sum(tmp > 0.5)/length(tmp)];
end
MedDNAmet = cat(1,MedDNAmet(1:4),MedDNAmet(5:8));
DNAmetscore = nanmean(MedDNAmet,1);

%%  average median DNAme
ymax = 0.75;
SizedFig(20,30);
hold on;
cond_names = {'no dox','5 days','irreversible','reactivated'};
xaxtmp = [0,1,2,3];
bar(xaxtmp,DNAmetscore(1:end),'FaceColor',[0.8,0.8,0.8]);
plot(xaxtmp,MedDNAmet(:,1:end),'o','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','w','Color',[0.9,0.2,0.2]);

set(gca,'XTick',[0,5,7,10,15],'XColor','k','YColor','k');

ylabel('median DNAme fraction');
set(gca,'XTick',xaxtmp(1:end),'XTickLabel',cond_names);
xtickangle(-45);
ylim([0,1]);