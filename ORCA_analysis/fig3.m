
%% directory for figure
figDIR = 'Z:\Taihei\project_ORCA\singleclones_ORCA\Analysis\fig1\';
%% re-aggregate data



dataselected = {'KRAB-nodox','KRAB-1day','HDAC4-nodox','HDAC4-1day','VP64-nodox','VP64-1day'};
ctrlidx = [1,1,3,3,5,5];

dfig3 = [];
dfig3.dmat = [];
dfig3.dmatNoIntp = [];
dfig3.dmatfilt = [];
dfig3.coordfilt = [];
dfig3.coordNoIntp = [];
dfig3.Rg = [];
dfig3.Rgfilt = [];
dfig3.rep = [];
dfig3.rep_dmatfilt = [];
dfig3.ctrlidx = [];
dfig3.ctrlidx_dmatfilt = [];
dfig3.dname = {};
dfig3.dname_dmatfilt = {};
for kk = 1:length(Objs)
   df1 = Objs(kk).df;
    drange = 1:length(df1);
    for ii = 1:length(drange)
        if sum(strcmp(df1(ii).dname,dataselected))
            dnum = size(df1(ii).dmatinterp,3);
            dfig3.dmat = cat(3,dfig3.dmat,df1(ii).dmatinterp);
            tmpRg = zeros(13,13,dnum);
            for jj = 4:15
                for mm = (jj+1):16
                    idxrange = jj:mm;
                    tmp = df1(ii).coordinterp(idxrange,:,:);
                    tmp = tmp - repmat(mean(tmp,1),[size(tmp,1), 1, 1]);
                    tmpRg(jj-3,mm-3,:) = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
                end
            end
            dfig3.Rg = cat(3,dfig3.Rg,tmpRg);
            dfig3.dmatNoIntp = cat(3,dfig3.dmatNoIntp,df1(ii).dmatfilt50p);
            dfig3.coordNoIntp = cat(3,dfig3.coordNoIntp,df1(ii).coordfilt50p);
            dfig3.dmatfilt = cat(3,dfig3.dmatfilt,df1(ii).dmatfilt);
            dfig3.coordfilt = cat(3,dfig3.coordfilt,df1(ii).coordfilt);
            for jj = 1:dnum
                dfig3.rep = cat(1,dfig3.rep,kk);            
                dfig3.ctrlidx = cat(1,dfig3.ctrlidx,ctrlidx(strcmp(df1(ii).dname,dataselected)));            
                dfig3.dname = [dfig3.dname;df1(ii).dname];            
            end
            
            tmpRg = zeros(13,13,size(df1(ii).dmatfilt,3));
            for jj = 4:15
                for mm = (jj+1):16
                    idxrange = jj:mm;
                    tmp = df1(ii).coordfilt(idxrange,:,:);
                    tmp = tmp - repmat(mean(tmp,1),[size(tmp,1), 1, 1]);
                    tmpRg(jj-3,mm-3,:) = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
                end
            end
            dfig3.Rgfilt = cat(3,dfig3.Rgfilt,tmpRg);
            for jj = 1:size(df1(ii).dmatfilt,3)
                dfig3.dname_dmatfilt = [dfig3.dname_dmatfilt;df1(ii).dname];            
                dfig3.rep_dmatfilt = cat(1,dfig3.rep_dmatfilt,kk);            
            end
        end
    end
end


%% -----------------------------------------------------------------
%% --------- calculate median distance map for each bioreps --------
%% -----------------------------------------------------------------
%% re-group data by experimental condition
% SizedFig(40,10);
diagidx = 4:16;
% diagidx = 1:13;
% diagidx = 1:19;
xlabels = {'–30','–20','–10','0','10','20','30'};
medmats = [];
medmats.dmat = [];
medmats.dmatNoIntp = [];
for ii = 1:length(dataselected)
%     SizedFig(60,10);
%     subplot(1,length(dataselected),ii);
    medmat = [];
    medmat2 = [];
    medmat3 = [];
    medmat4 = [];
    for jj = 1:length(Objs)
        subplot(1,length(Objs),jj);
        tmp = dfig3.dmat(:,:,strcmp(dfig3.dname,dataselected(ii)) & dfig3.rep == jj);
        tmp = nanmedian(tmp,3);
        medmat = cat(3,medmat,tmp);
%         if ~isnan(mean(tmp(:)))
%             imagetriu(tmp(diagidx,diagidx),200,300,flipud(jet),0);
%             title(dataselected(ii));
%             caxis([200,300]);
%             cbar = colorbar;
%             cbar.Label.String = 'nm';
%             cbar.Label.FontSize = 15;
%             set(gca,'XTick',1:2:13,'XTickLabels',xlabels,'FontSize',10);
%         end
        tmp = dfig3.dmatfilt(:,:,strcmp(dfig3.dname_dmatfilt,dataselected(ii)) & dfig3.rep_dmatfilt == jj);
        tmp = nanmedian(tmp,3);
        medmat2 = cat(3,medmat2,tmp);
        tmp = dfig3.Rg(:,:,strcmp(dfig3.dname,dataselected(ii)) & dfig3.rep == jj);
        tmp = nanmedian(tmp,3);
        medmat3 = cat(3,medmat3,tmp);
        tmp = dfig3.Rgfilt(:,:,strcmp(dfig3.dname_dmatfilt,dataselected(ii)) & dfig3.rep_dmatfilt == jj);
        tmp = nanmedian(tmp,3);
        medmat4 = cat(3,medmat4,tmp);
    end
    medmats(ii).dmat = medmat;
    medmats(ii).dmatfilt = medmat2(1:end-1,1:end-1,:);
    medmats(ii).Rg = medmat3;
    medmats(ii).Rgfilt = medmat4;
end
%% subtracted median distancemap 
% SizedFig(40,10);
diagidx = 4:16;
% diagidx = 1:13;
% diagidx = 1:19;
xlabels = {'–30','–20','–10','0','10','20','30'};

for ii = 1:length(dataselected)
    SizedFig(40,10);
    subplot(1,length(dataselected),ii);
    for jj = 1:5
        subplot(1,5,jj);
        tmp = dfig3.dmat(:,:,strcmp(dfig3.dname,dataselected(ii)) & dfig3.rep == jj);
        tmp = nanmedian(tmp,3);
        tmp2 = dfig3.dmat(:,:,strcmp(dfig3.dname,dataselected(ctrlidx(ii))) & dfig3.rep == jj);
        tmp2 = nanmedian(tmp2,3);
%         tmp = dfig3.dmatfilt(:,:,strcmp(dfig3.dname_dmatfilt,dataselected(ii)) & dfig3.rep_dmatfilt == jj);
%         tmp = nanmedian(tmp,3);
%         tmp2 = dfig3.dmatfilt(:,:,strcmp(dfig3.dname_dmatfilt,dataselected(ctrlidx(ii))) & dfig3.rep_dmatfilt == jj);
%         tmp2 = nanmedian(tmp2,3);
        tmp = tmp - tmp2;
        if ~isnan(mean(tmp(:)))
            imagetriu(tmp(diagidx,diagidx),-100,100,flipud(bluewhiteredw0),0);
            title(dataselected(ii));
            caxis([-100,100]);
            cbar = colorbar;
            cbar.Label.String = 'nm';
            cbar.Label.FontSize = 15;
            set(gca,'XTick',1:2:13,'XTickLabels',xlabels,'FontSize',10);
        end
    end
end
%% AVG -- median distancemap 
% SizedFig(60,10);
fh = SizedFig(100,20);

diagidx = 1:19;
% diagidx = 4:16;
% diagidx = 1:13;
% diagidx = [1:3,4,7,10,13,16,17:19];

vizidx = [2,4,6];%1:length(medmats);%
% vizidx = [9,10,11];%1:length(medmats);%
% xlabels = {'–30','–20','–10','0','10','20','30'};
xlabels = {'–75','–30','–15','0','15','30','75'};
xlabelspos = [1,4,7,10,13,16,19];
for ii = 1:length(vizidx)
    subplot(1,length(vizidx),ii);
    medmat = nanmean(medmats(vizidx(ii)).dmatfilt,3);
    imagetriu(rot90(medmat(diagidx,diagidx),2),230,330,flipud(jet),0);
%     medmat = nanmean(medmats(vizidx(ii)).Rg,3);
%     imagetriu(medmat(diagidx,diagidx),50,180,flipud(jet),0);
    title(dataselected(vizidx(ii)));
    caxis([230,330]);
    cbar = colorbar;
    cbar.Label.String = 'nm';
    cbar.Label.FontSize = 15;
    set(gca,'XTick',xlabelspos,'XTickLabels',xlabels,'FontSize',10,...
        'XColor','k','YColor','k');
    xlabel('kb');
%     saveas(fh,[figDIR,'avg_medmat_',dataselected{vizidx(ii)},'.pdf']);
end

%% AVG -- subtraction of median distancemap xreversed
% SizedFig(60,10);
fh = SizedFig(100,20);

% diagidx = [1,2,3,4,7,10,13,16,17,18,19];%1:19;
% diagidx = 1:13;
diagidx = 4:16;
diagidx = 1:19;

vizidx = [2,4,6];
% vizidx = [7,8,10,11,13];
% vizidx = 1:length(medmats);%[1,3,4,5];
% vizidx = [2,3, 7,8];
% vizidx = [2,3, 7,8];
% vizidx = [9,10,11];
% xlabels = {'–75','–30','–15','0','15','30','75'};
% xlabelspos = [1,4,7,10,13,16,19];
xlabels = {'–30','–15','0','15','30'};
xlabelspos = [1,4,7,10,13];
for ii = 1:length(vizidx)
    subplot(1,length(vizidx),ii);
%     medmat = nanmean(medmats(vizidx(ii)).dmat - medmats(ctrlidx(vizidx(ii))).dmat,3);
    medmat = nanmean(medmats(vizidx(ii)).dmatfilt - medmats(ctrlidx(vizidx(ii))).dmatfilt,3);
    imagetriu(rot90(medmat(diagidx,diagidx),2),-100,100,flipud(bluewhiteredw0),0);
%     medmat = nanmean(medmats(vizidx(ii)).Rgfilt - medmats(ctrlidx(vizidx(ii))).Rgfilt,3);
%     imagetriu(medmat(diagidx,diagidx),-50,50,flipud(bluewhiteredw0),0);
    title(dataselected(vizidx(ii)));
    caxis([-100,100]);
    cbar = colorbar;
    cbar.Label.String = 'nm';
    cbar.Label.FontSize = 15;
    set(gca,'XTick',xlabelspos,'XTickLabels',xlabels,'FontSize',10);
    xlabel('kb');
%     saveas(fh,[figDIR,'subtraction_medmat_',dataselected{vizidx(ii)},'.pdf']);
end

%% local compaction by Rg Ratio; bar mean, xreversed
{'KRAB-nodox','KRAB-1day','HDAC4-nodox','HDAC4-1day','VP64-nodox','VP64-1day'};
ii=2;
diagidx = 1:13;
medmatd = medmats(ii).Rg(diagidx,diagidx,:);
medmatc = medmats(ctrlidx(ii)).Rg(diagidx,diagidx,:);
localRgD = [];
localRgC = [];
for jj = 1:(size(medmatd,1)-4)
    localRgD(jj,:) = medmatd(jj,jj+4,:);
    localRgC(jj,:) = medmatc(jj,jj+4,:);
end
Rgratio = flipud(localRgD./localRgC);
avgRgRatio = nanmean(Rgratio,2);
stdRgRatio = nanstd(Rgratio,[],2) / sqrt(sum(~isnan(Rgratio(1,:))));

SizedFig(30,25);
xax = linspace(-floor(size(localRgD,1)/2)*5,floor(size(localRgD,1)/2)*5,size(localRgD,1));
lh = plot(repmat(xax,[2,1]),[ones(size(localRgD,1),1),avgRgRatio]','b-','linewidth',25);
for ii = 1:length(lh)
    lh(ii).Color=[0,0,1,0.5];
end
hold on;
plot(repmat(xax,[2,1]),[avgRgRatio-stdRgRatio,avgRgRatio+stdRgRatio]','k-','linewidth',1);
% for ii = 1:7
%     plot(xax+(ii-5.5)/5,params(:,ii),'bo','MarkerFaceColor','w','MarkerSize',6);
% end
plot([-30,30],[1.0,1.0],'k-');
plot([0,0],[0.2,2],'k--');
ylim([0.86,1.2]); % for others
% ylim([0.82,1.26]); % for VP64
xlim([-30,30]);
xlabel('position [kb]');
ylabel('Rg ratio (treat/nodox)');
set(gca,'FontSize',15,'XColor','k','YColor','k');
box on;


%% histogram of radius of gyration from specified biorep
for ii = 1:length(Objs)
    for kk = 1:length(Objs(ii).df)
        Objs(ii).df(kk).rep = ii;
    end
end
for ii = 1:length(Objs)
    for kk = 1:length(Objs(ii).df)
        dftmp = Objs(ii).df(kk);
        tmpRg = zeros(13,13,size(dftmp.coordinterp,3));
        for jj = 4:15
            for mm = (jj+1):16
                idxrange = jj:mm;
                tmp = dftmp.coordinterp(idxrange,:,:);
                tmp = tmp - repmat(mean(tmp,1),[size(tmp,1), 1, 1]);
                tmpRg(jj-3,mm-3,:) = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
            end
        end
        Objs(ii).df(kk).Rg = tmpRg;
    end
end
dataselected = {'KRAB-nodox','KRAB-1day','HDAC4-nodox','HDAC4-1day','VP64-nodox','VP64-1day'};
dfall = [];
for idxobj = 1:length(Objs)
    for jj = 1:length(Objs(idxobj).df)
        if sum(strcmp(Objs(idxobj).df(jj).dname,dataselected))
            dfall = cat(2,dfall,Objs(idxobj).df(jj));
        end
    end
end
for ii = 1:length(dfall)
    allpairs = [];
    tmpidx = ~isnan(reshape(dfall(ii).coordfilt50p(10,1,:),[],1));
    tmp = dfall(ii).dmatinterp(4:16,4:16,tmpidx);
    tmp_triu = ~triu(ones(size(tmp,1),size(tmp,2)))';
    disp([dfall(ii).dname,':',num2str(size(tmp,3))]);
    for i = 1:size(tmp,3)
        tmpmat = tmp(:,:,i);
        allpairs = cat(1,allpairs,tmpmat(tmp_triu)');
    end
    dfall(ii).allpairs = allpairs;
    idx = 1:size(tmp,3);%randsample(size(tmp,3),267);
    dfall(ii).pairs = allpairs(idx,:);
end
for ii = 1:length(dfall)
    disp([num2str(ii),':',dfall(ii).dname]);
end
%%
cmap = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980];
SizedFig(15,15);
hold on;
idx1 = 3;%1;% % krab 1 day
idx2 = 4;%2;% % krab 1 day
%  idx1 = 7;%13;% % irr
%  idx2 = 8;%14;% % irr
%  idx1 = 16;%18;%
%  idx2 = 17;%19;%
binnum = 20;
binr = [50,400];
tmp1 = reshape(dfall(idx1).Rg(1,13,:),[],1);
h1 = histogram(tmp1,binnum,'BinLimits',binr,'normalization','probability',...
    'EdgeAlpha',0.0,'FaceColor',cmap(1,:));
stairs([h1.BinEdges,h1.BinEdges(end)],[h1.Values,h1.Values(end),0],'k-');

tmp2 = reshape(dfall(idx2).Rg(1,13,:),[],1);
h2 = histogram(tmp2,binnum,'BinLimits',binr,'normalization','probability',...
    'EdgeAlpha',0.0,'FaceColor',cmap(2,:));
stairs([h2.BinEdges,h2.BinEdges(end)],[h2.Values,h2.Values(end),0],'k-');

plot([median(tmp1),median(tmp1)],[0,0.12],'--','Color',cmap(1,:));
plot([median(tmp2),median(tmp2)],[0,0.12],'--','Color',cmap(2,:));
ylim([0,0.2])
set(gca,'XColor','k','YColor','k');
median(tmp1)
median(tmp2)
ranksum(tmp1,tmp2)
box on;
xlabel('Radius of gyration');
ylabel('probability');




%%
DIR = '/Users/tfuji/Documents/tempanalysis/_OnlyReporterPos/fig3/submat/';
vizidx = [6];
ii = 1;
diagidx = 4:16;
medmat = nanmean(medmats(vizidx(ii)).dmat - medmats(ctrlidx(vizidx(ii))).dmat,3);
writematrix(medmat(diagidx,diagidx),[DIR,'dmatdiff_',dataselected{vizidx(ii)},'.csv']);
%% calcaverage neighbor length
neighborlength = [];
for ii = 1:length(medmats)
    tmp = medmats(ii).dmat(diagidx,diagidx,:);
    for jj = 1:size(tmp,3)
        medmat = tmp(:,:,jj);
        neighborlength = [neighborlength,sum(sum((triu(medmat,1)-triu(medmat,2))))/12];
    end
end
neighborlength = neighborlength(~isnan(neighborlength));
mean(neighborlength)
%%
figure;
imagesc(medmat(diagidx,diagidx));
axis image;
title(dataselected(vizidx(ii)));
caxis([-30,30]);
cbar = colorbar;
cbar.Label.String = 'nm';
cbar.Label.FontSize = 15;
set(gca,'XTick',1:2:13,'XTickLabels',xlabels,'FontSize',10);
xlabel('kb');


%%
figure; 
hold on; 
pval = [1.0370378512793816,0.8667000040178601,1.010993160477617,0.9054793188036305,0.8871767809532564,0.7804728499237115,0.8992139386739062,0.8579843088855152,1.0052844186119696,0.8487996568417516,1.044876372408926,1.0137938976928014];
plot(pval,'bo'); 
pval = [0.8103175920744072,1.0510810453012258,1.0339708802650125,0.843263403566806,0.9395586846228271,0.7327399132915315,0.9971397761090393,0.8676781870661499,0.9507136084060674,0.9811129275221269,0.9655082485099366,1.0329214828303517];
plot(pval,'bo'); 
plot([1,12],[1,1],'k--');
plot([6.5,6.5],[0.6,1.1],'k--');
ylim([0.6,1.1]);
%% adjacent segment distance; subtraction
SizedFig(15,20);
hold on;
xax = [-75,-60,-45,-30:5:30,45,60,75]'; 
xax = (-30:5:30)';
diagidx = 4:16;
% diagidx = 1:13;
xax = (xax(1:end-1)+xax(2:end))/2;
xlabels = {'–30','–20','–10','0','10','20','30'};

visidx = [4,3,5];%reactivate,KRAB5days,irrversible
visidx = [2,3];%krab 1day,5days
% visidx = [7,8];%krab147 1day,5days
% visidx = [3,8];%krab 5days, krab147 5days
visidx = [2,3,7,8];%1:length(dataselected);
% visidx = [4,3,5];%reactivate,KRAB5days,irrversible
cmap = parula(length(visidx)+2);
cmap = cmap(2:end-1,:);
% cmap = [0.6 0.1 0.6; 0.7 0.4 0.1;0.1 0.6 0.7];
plot([-80,80],[0,0],'k--','LineWidth',1);
for ii = 1:length(visidx)
    adjdist = nan(length(xax),5);
    medmat = medmats(visidx(ii)).dmat(diagidx,diagidx,:) - medmats(ctrlidx(visidx(ii))).dmat(diagidx,diagidx,:);
%     medmat = medmats(visidx(ii)).Rg(diagidx,diagidx,:) - medmats(ctrlidx(visidx(ii))).Rg(diagidx,diagidx,:);
%     medmat = medmats(visidx(ii)).dmatfilt - medmats(ctrlidx(visidx(ii))).dmatfilt;
    for jj = 1:(size(medmat,1)-1)
        adjdist(jj,:) = medmat(jj,jj+1,:);
    end
    avgdist = nanmean(adjdist,2);
    stddist = nanstd(adjdist,[],2) / sqrt(sum(~isnan(adjdist(1,:))));
    plot(xax,avgdist,'o-','Color',cmap(ii,:),'LineWidth',2);
    plot([xax,xax]',[avgdist-stddist,avgdist+stddist]','-','Color',cmap(ii,:),'LineWidth',1);
end
% xlim([-80,80]);
xlim([-32.5,32.5]);
% ylim([-20,20]);
set(gca,'FontSize',10);
xlabel('position [kb]');
ylabel('difference from nodox [nm]');

%% adjacent segment distance; absolute distance
SizedFig(15,20);
hold on;
xax = [-75,-60,-45,-30:5:30,45,60,75]';
perkb = repmat(xax(2:end) - xax(1:end-1),[1 5]);
xax = (-30:5:30)';
xax = (xax(1:end-1)+xax(2:end))/2;
% diagidx = 4:16;
diagidx = 4:16;
xlabels = {'–30','–20','–10','0','10','20','30'};

visidx = [2,3,4,5,7,8];%1:length(dataselected);
% visidx = [7,8,2,3];%1:length(dataselected);
visidx = [1,4,3,5];%1:length(dataselected);
% visidx = [1,5];%1:length(dataselected);
cmap = parula(length(visidx)+2);
cmap = cmap(2:end-1,:);
% cmap = [0.6 0.1 0.6; 0.7 0.4 0.1;0.1 0.6 0.7];
% plot([-80,80],[0,0],'k--','LineWidth',1);    
for ii = 1:length(visidx)
    adjdist = nan(length(xax),5);
    medmat = medmats(visidx(ii)).dmat(diagidx,diagidx,:);%;%dmatfilt(4:16,4:16,:);%
    for jj = 1:(size(medmat,1)-1)
        adjdist(jj,:) = medmat(jj,jj+1,:);%/perkb(jj);
    end
    avgdist = nanmean(adjdist,2);
    stddist = nanstd(adjdist,[],2) / sqrt(sum(~isnan(adjdist(1,:))));
    plot(xax,avgdist,'o-','Color',cmap(ii,:),'LineWidth',2);
    plot([xax,xax]',[avgdist-stddist,avgdist+stddist]','-','Color',cmap(ii,:),'LineWidth',1);
end
xlim([-80,80]);
% xlim([-32.5,32.5]);
% ylim([-45,25]);
% set(gca, 'YDir','reverse','FontSize',10);
xlabel('position [kb]');
ylabel('difference from nodox [nm]');

%% end-to-end distance
SizedFig(15,20);
diagidx = 1:19;
xax = [-75,-60,-45,-30:5:30,45,60,75]';
xax = xax(11:end)*2;
xlabels = {'–30','–20','–10','0','10','20','30'};

visidx = 1:length(dataselected);%[1,3,4,5];%
cmap = parula(length(visidx)+2);
cmap = cmap(2:end-1,:);
for ii = 1:length(visidx)
    adjdist = nan(length(xax),5);
    medmat = medmats(ctrlidx(visidx(ii))).dmat - medmats(visidx(ii)).dmat;
%     medmat = medmats(visidx(ii)).dmat;
    center = round(size(medmat,1)/2);
    for jj = 1:(center-1)
        adjdist(jj,:) = medmat(center-jj,center+jj,:);
    end
    avgdist = nanmean(adjdist,2);
    stddist = nanstd(adjdist,[],2) / sqrt(sum(~isnan(adjdist(1,:))));
    plot(xax,avgdist,'o-','Color',cmap(ii,:),'LineWidth',2);
%     plot([xax,xax]',[avgdist-stddist,avgdist+stddist]','-','Color',cmap(ii,:),'LineWidth',1);
    hold on;
end
legend(dataselected(visidx),'location','eastoutside');
plot([0,160],[0,0],'k--','LineWidth',1);    
xlim([0,160]);    
xlabel('end-to-end distance [kb]');
ylabel('nodox - dox [nm]');
    
%% distance to reporter region
SizedFig(20,20);
hold on;
diagidx = 1:19;
xax = [-75,-60,-45,-30:5:30,45,60,75]';
xlabels = {'–30','–20','–10','0','10','20','30'};

visidx = [2,3,4,5];%1:length(dataselected);
cmap = parula(length(visidx)+2);
cmap = cmap(2:end-1,:);
plot([-80,80],[0,0],'k--','LineWidth',1);    
for ii = 1:length(visidx)
    adjdist = nan(length(xax),5);
    medmat = medmats(ctrlidx(visidx(ii))).dmat - medmats(visidx(ii)).dmat;
    for jj = 1:size(medmat,1)
        adjdist(jj,:) = medmat(jj,10,:);
    end
    avgdist = nanmean(adjdist,2);
    stddist = nanstd(adjdist,[],2) / sqrt(sum(~isnan(adjdist(1,:))));
    plot(xax,avgdist,'o-','Color',cmap(ii,:),'LineWidth',2);
    plot([xax,xax]',[avgdist-stddist,avgdist+stddist]','-','Color',cmap(ii,:),'LineWidth',1);
end
xlim([-80,80]);
xlabel('position [kb]');
ylabel('nodox - dox [nm]');
        
%%

avgRgFC = [];
stdRgFC = [];
visidx = 1:length(dataselected);
for ii = 1:length(visidx)
    medmat = medmats(visidx(ii)).Rg(1,13,:) ./ medmats(ctrlidx(visidx(ii))).Rg(1,13,:);
    avgRgFC = [avgRgFC;nanmean(medmat,3)];
    stdRgFC = [stdRgFC;nanstd(medmat,[],3)/sqrt(sum(~isnan(medmat)))];
end




%% KRAB 5days closest to medmat
diagidx = 4:16;
dmat1 = cat(3,Objs(1).df(3).dmatinterp(diagidx,diagidx,:),Objs(2).df(3).dmatinterp(diagidx,diagidx,:),...
    Objs(3).df(3).dmatinterp(diagidx,diagidx,:),Objs(5).df(3).dmatinterp(diagidx,diagidx,:));
tmpcoord = cat(3,Objs(1).df(3).coordinterp(diagidx,:,:),Objs(2).df(3).coordinterp(diagidx,:,:),...
    Objs(3).df(3).coordinterp(diagidx,:,:),Objs(5).df(3).coordinterp(diagidx,:,:));
refmat = nanmean(medmats(3).dmat(diagidx,diagidx,:),3);
difffrommed = reshape(sum(sum((dmat1 - refmat).^2,1),2),[],1);
% difffrommed = reshape(sum(sum(abs(dmat1 - refmat),1),2),[],1);
idxclosest = find(difffrommed == min(difffrommed))
%
figure; 
subplot(1,2,1);
imagetriu(refmat,200,300,flipud(jet));
subplot(1,2,2);
imagetriu(dmat1(:,:,idxclosest),0,300,flipud(jet));
az = -161;
el = -12;
xangle = pi/4; yangle = -pi/8; zangle = pi;
az = 31;
el = 19;
repcoordkrab = tmpcoord(:,:,idxclosest);
%% KRAB no dox closest to medmat
diagidx = 4:16;
dmat1 = cat(3,Objs(1).df(1).dmatinterp(diagidx,diagidx,:),Objs(2).df(1).dmatinterp(diagidx,diagidx,:),...
    Objs(3).df(1).dmatinterp(diagidx,diagidx,:),Objs(4).df(1).dmatinterp(diagidx,diagidx,:),Objs(5).df(1).dmatinterp(diagidx,diagidx,:));
tmpcoord = cat(3,Objs(1).df(1).coordinterp(diagidx,:,:),Objs(2).df(1).coordinterp(diagidx,:,:),...
    Objs(3).df(1).coordinterp(diagidx,:,:),Objs(4).df(1).coordinterp(diagidx,:,:),Objs(5).df(1).coordinterp(diagidx,:,:));
refmat = nanmean(medmats(1).dmat(diagidx,diagidx,:),3);
difffrommed = reshape(sum(sum((dmat1 - refmat).^2,1),2),[],1);
% difffrommed = reshape(sum(sum(abs(dmat1 - refmat),1),2),[],1);
idxclosest = find(difffrommed == min(difffrommed))
%
figure; 
subplot(1,2,1);
imagetriu(refmat,200,300,flipud(jet));
subplot(1,2,2);
imagetriu(dmat1(:,:,idxclosest),0,300,flipud(jet));
az = 31;
el = 19;
xangle = 0; yangle = 0; zangle = 0;
repcoordnodox = tmpcoord(:,:,idxclosest);

%% visualize the 3D structure
cntr = (idxrange(end) - idxrange(1))/2 + 1;
cmap = turbo(size(tmpcoord,1)+2);
cmap = cmap(2:end-1,:);
boxhw = 200;
SizedFig(30,40); 
tmp = repcoordnodox; 
tmp = tmp(~isnan(tmp(:,1)),:); 
polycenter = mean(tmp,1);
tmp = tmp - polycenter;
tmp = (rotz(zangle)*roty(yangle)*rotx(xangle)*tmp')';
title('+ dox');    
minmax = [min(tmp,[],1) - 0.3*abs(max(tmp,[],1) - min(tmp,[],1));...
    max(tmp,[],1) + 0.3*abs(max(tmp,[],1) - min(tmp,[],1))];

if sum(abs(minmax(1,:) - minmax(2,:))) ~= 0 && ~isnan(sum(abs(minmax(1,:) - minmax(2,:))))
    tmpinterpx = [];
    tmpinterpy = [];
    tmpinterpz = [];
    tmpinterpx = cat(1,tmpinterpx,...
        interp1(1:size(tmp,1),tmp(:,1),1:0.1:size(tmp,1),'linear')');
    tmpinterpy = cat(1,tmpinterpy,...
        interp1(1:size(tmp,1),tmp(:,2),1:0.1:size(tmp,1),'linear')');
    tmpinterpz = cat(1,tmpinterpz,...
        interp1(1:size(tmp,1),tmp(:,3),1:0.1:size(tmp,1),'linear')');
    xrange = [max(tmpinterpx),min(tmpinterpx)];
    yrange = [max(tmpinterpy),min(tmpinterpy)];
    zrange = [max(tmpinterpz),min(tmpinterpz)];
    widthMax = max([xrange(1)-xrange(2),yrange(1)-yrange(2),zrange(1)-zrange(2)])*1.15/2;
    xbound = [sum(xrange)/2-widthMax,sum(xrange)/2+widthMax];
    ybound = [sum(yrange)/2-widthMax,sum(yrange)/2+widthMax];
    zbound = [sum(zrange)/2-widthMax,sum(zrange)/2+widthMax];
    radius = 100*widthMax/1000;
    xbound = [-boxhw,boxhw];
    ybound = [-boxhw,boxhw];
    zbound = [-boxhw,boxhw];
    radius = 20;
    clf;%
    hold on;
    [xsphere,ysphere,zsphere] = sphere; 
    xsphere = xsphere*radius; ysphere = ysphere*radius; zsphere = zsphere*radius;
    for i =1:length(tmp)
        % --- plot spheres --- %
        surf(xsphere+tmp(i,1),ysphere+tmp(i,2),zsphere+tmp(i,3),'EdgeAlpha',0,'FaceColor',cmap(i,:));
        if i ~= length(tmp)
            % --- plot spline interpolation --- %
            pltrange = (1+10*(i-1)):(1+10*i);
            h = plot3t(tmpinterpx(pltrange),tmpinterpy(pltrange),tmpinterpz(pltrange),...
                30*widthMax/1000,cmap(i,:));
            % --- optimize axis --- %
            set(h, 'FaceLighting','phong','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
            set(h,'EdgeAlpha',0);
            material shiny; 
        end
    end

    % --- add texts --- %
%         texts = {'  pEF-mCit'};
%         txtidx = cntr;
%         for i = 1:length(texts)
%             text(tmp(txtidx,1),tmp(txtidx,2),tmp(txtidx,3),texts{i},'FontSize',15)
%         end
    title(['spot#',num2str(ii)]);


    % --- optimize axis --- %
    axis equal;
    grid on
    set(gca,'Projection','perspective','Box','on','BoxStyle','full',...
        'FontSize',10)
    xlim(xbound); ylim(ybound); zlim(zbound);
    xlabel('X [nm]');
    ylabel('Y [nm]');
    zlabel('Z [nm]');
    set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
%         view(az(21),el);
    view(az,el);
%         view([7.9927 28.5909]);
    camlight('right');
    cntr = 0;
    pause(1/10);
end
box off;
%% visualize the 3D structure
cntr = (idxrange(end) - idxrange(1))/2 + 1;

tmp = repcoordnodox; 
xangle = 0; yangle = 0; zangle = 0;
% tmp = repcoordkrab; 
% xangle = pi/4; yangle = -pi/8; zangle = pi;

cmap = turbo(size(tmp,1)*2+2);
cmap = cmap(2:end-1,:);
% cmap(1:2:end,:) = repmat([0.5,0.5,0.5],[13,1]);
boxhw = 200;

SizedFig(30,40); 
tmp = tmp(~isnan(tmp(:,1)),:); 
polycenter = mean(tmp,1);
tmp = tmp - polycenter;
tmp = (rotz(zangle)*roty(yangle)*rotx(xangle)*tmp')';
title('+ dox');    
minmax = [min(tmp,[],1) - 0.3*abs(max(tmp,[],1) - min(tmp,[],1));...
    max(tmp,[],1) + 0.3*abs(max(tmp,[],1) - min(tmp,[],1))];

if sum(abs(minmax(1,:) - minmax(2,:))) ~= 0 && ~isnan(sum(abs(minmax(1,:) - minmax(2,:))))
    tmpinterpx = [];
    tmpinterpy = [];
    tmpinterpz = [];
    tmpinterpx = cat(1,tmpinterpx,...
        interp1(1:size(tmp,1),tmp(:,1),1:0.1:size(tmp,1),'linear')');
    tmpinterpy = cat(1,tmpinterpy,...
        interp1(1:size(tmp,1),tmp(:,2),1:0.1:size(tmp,1),'linear')');
    tmpinterpz = cat(1,tmpinterpz,...
        interp1(1:size(tmp,1),tmp(:,3),1:0.1:size(tmp,1),'linear')');
    xrange = [max(tmpinterpx),min(tmpinterpx)];
    yrange = [max(tmpinterpy),min(tmpinterpy)];
    zrange = [max(tmpinterpz),min(tmpinterpz)];
    widthMax = max([xrange(1)-xrange(2),yrange(1)-yrange(2),zrange(1)-zrange(2)])*1.15/2;
    xbound = [sum(xrange)/2-widthMax,sum(xrange)/2+widthMax];
    ybound = [sum(yrange)/2-widthMax,sum(yrange)/2+widthMax];
    zbound = [sum(zrange)/2-widthMax,sum(zrange)/2+widthMax];
    radius = 100*widthMax/1000;
    xbound = [-boxhw,boxhw];
    ybound = [-boxhw,boxhw];
    zbound = [-boxhw,boxhw];
    radius = 20;
    clf;%
    hold on;
    [xsphere,ysphere,zsphere] = sphere; 
    xsphere = xsphere*radius; ysphere = ysphere*radius; zsphere = zsphere*radius;
    for i =1:length(tmp)
        % --- plot spheres --- %
        surf(xsphere+tmp(i,1),ysphere+tmp(i,2),zsphere+tmp(i,3),'EdgeAlpha',0,'FaceColor',cmap(i*2-1,:));
        if i ~= length(tmp)
            % --- plot spline interpolation --- %
            pltrange = (1+10*(i-1)):(1+10*i);
            h = plot3t(tmpinterpx(pltrange),tmpinterpy(pltrange),tmpinterpz(pltrange),...
                30*widthMax/1000,cmap(i*2,:));
            % --- optimize axis --- %
            set(h, 'FaceLighting','phong','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
            set(h,'EdgeAlpha',0);
            material shiny; 
        end
    end

    % --- add texts --- %
%         texts = {'  pEF-mCit'};
%         txtidx = cntr;
%         for i = 1:length(texts)
%             text(tmp(txtidx,1),tmp(txtidx,2),tmp(txtidx,3),texts{i},'FontSize',15)
%         end

    % --- optimize axis --- %
    axis equal;
    grid on
    set(gca,'Projection','perspective','Box','off','BoxStyle','full',...
        'FontSize',10)
    xlim(xbound); ylim(ybound); zlim(zbound);
    xlabel('X [nm]');
    ylabel('Y [nm]');
    zlabel('Z [nm]');
    set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
%         view(az(21),el);
    view(az,el);
%         view([7.9927 28.5909]);
    camlight('right');
    cntr = 0;
    pause(1/10);
end
box off;
%%

rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)];
roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)];
rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1]; 

