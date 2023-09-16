
%% directory for figure
figDIR = 'Z:\Taihei\project_ORCA\singleclones_ORCA\Analysis\fig1\';
%% re-aggregate data

dataselected =  {'ctrl','difday1','difday2','difday3'};
ctrlidx = [1,1,1,1];

dfig5 = [];
dfig5.dmat = [];
dfig5.dmatNoIntp = [];
dfig5.dmatfilt = [];
dfig5.coordfilt = [];
dfig5.coordNoIntp = [];
dfig5.Rg = [];
dfig5.Rgfilt = [];
dfig5.rep = [];
dfig5.rep_dmatfilt = [];
dfig5.ctrlidx = [];
dfig5.ctrlidx_dmatfilt = [];
dfig5.dname = {};
dfig5.dname_dmatfilt = {};
for kk = 1:length(Objs)
   df1 = Objs(kk).df;
    drange = 1:length(df1);
    for ii = 1:length(drange)
        if sum(strcmp(df1(ii).dname,dataselected))
            dnum = size(df1(ii).dmatinterp,3);
            dfig5.dmat = cat(3,dfig5.dmat,df1(ii).dmatinterp);
            tmpRg = zeros(13,13,dnum);
            for jj = 4:15
                for mm = (jj+1):16
                    idxrange = jj:mm;
                    tmp = df1(ii).coordinterp(idxrange,:,:);
                    tmp = tmp - repmat(mean(tmp,1),[size(tmp,1), 1, 1]);
                    tmpRg(jj-3,mm-3,:) = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
                end
            end
            dfig5.Rg = cat(3,dfig5.Rg,tmpRg);
            dfig5.dmatNoIntp = cat(3,dfig5.dmatNoIntp,df1(ii).dmatfilt50p);
            dfig5.coordNoIntp = cat(3,dfig5.coordNoIntp,df1(ii).coordfilt50p);
            dfig5.dmatfilt = cat(3,dfig5.dmatfilt,df1(ii).dmatfilt);
            dfig5.coordfilt = cat(3,dfig5.coordfilt,df1(ii).coordfilt);
            for jj = 1:dnum
                dfig5.rep = cat(1,dfig5.rep,kk);            
                dfig5.ctrlidx = cat(1,dfig5.ctrlidx,ctrlidx(strcmp(df1(ii).dname,dataselected)));            
                dfig5.dname = [dfig5.dname;df1(ii).dname];            
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
            dfig5.Rgfilt = cat(3,dfig5.Rgfilt,tmpRg);
            for jj = 1:size(df1(ii).dmatfilt,3)
                dfig5.dname_dmatfilt = [dfig5.dname_dmatfilt;df1(ii).dname];            
                dfig5.rep_dmatfilt = cat(1,dfig5.rep_dmatfilt,kk);            
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
    SizedFig(60,10);
    subplot(1,length(dataselected),ii);
    medmat = [];
    medmat2 = [];
    medmat3 = [];
    medmat4 = [];
    for jj = 1:length(Objs)
        subplot(1,length(Objs),jj);
        tmp = dfig5.dmat(:,:,strcmp(dfig5.dname,dataselected(ii)) & dfig5.rep == jj);
        tmp = nanmedian(tmp,3);
        medmat = cat(3,medmat,tmp);
        if ~isnan(mean(tmp(:)))
            imagetriu(tmp(diagidx,diagidx),200,300,flipud(jet),0);
            title(dataselected(ii));
            caxis([200,300]);
            cbar = colorbar;
            cbar.Label.String = 'nm';
            cbar.Label.FontSize = 15;
            set(gca,'XTick',1:2:13,'XTickLabels',xlabels,'FontSize',10);
        end
        tmp = dfig5.dmatfilt(:,:,strcmp(dfig5.dname_dmatfilt,dataselected(ii)) & dfig5.rep_dmatfilt == jj);
        tmp = nanmedian(tmp,3);
        medmat2 = cat(3,medmat2,tmp);
        tmp = dfig5.Rg(:,:,strcmp(dfig5.dname,dataselected(ii)) & dfig5.rep == jj);
        tmp = nanmedian(tmp,3);
        medmat3 = cat(3,medmat3,tmp);
        tmp = dfig5.Rgfilt(:,:,strcmp(dfig5.dname_dmatfilt,dataselected(ii)) & dfig5.rep_dmatfilt == jj);
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
diagidx = 1:13;
diagidx = 1:19;
xlabels = {'–30','–20','–10','0','10','20','30'};

for ii = 1:length(dataselected)
    SizedFig(40,10);
    subplot(1,length(dataselected),ii);
    for jj = 1:2
        subplot(1,2,jj);
        tmp = dfig5.dmat(:,:,strcmp(dfig5.dname,dataselected(ii)) & dfig5.rep == jj);
        tmp = nanmedian(tmp,3);
        tmp2 = dfig5.dmat(:,:,strcmp(dfig5.dname,dataselected(ctrlidx(ii))) & dfig5.rep == jj);
        tmp2 = nanmedian(tmp2,3);
        tmp = dfig5.dmatfilt(:,:,strcmp(dfig5.dname_dmatfilt,dataselected(ii)) & dfig5.rep_dmatfilt == jj);
        tmp = nanmedian(tmp,3);
        tmp2 = dfig5.dmatfilt(:,:,strcmp(dfig5.dname_dmatfilt,dataselected(ctrlidx(ii))) & dfig5.rep_dmatfilt == jj);
        tmp2 = nanmedian(tmp2,3);
        tmp = tmp - tmp2;
        if ~isnan(mean(tmp(:)))
            imagetriu(tmp(diagidx,diagidx),-50,50,flipud(bluewhiteredw0),0);
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

vizidx = [1,2,3,4];%1:length(medmats);%
% vizidx = [9,10,11];%1:length(medmats);%
% xlabels = {'–30','–20','–10','0','10','20','30'};
xlabels = {'–75','–30','–15','0','15','30','75'};
xlabelspos = [1,4,7,10,13,16,19];
for ii = 1:length(vizidx)
    subplot(1,length(vizidx),ii);
    medmat = nanmean(medmats(vizidx(ii)).dmatfilt,3);
    imagetriu(medmat(diagidx,diagidx),150,250,flipud(jet),0);
%     medmat = nanmean(medmats(vizidx(ii)).Rg,3);
%     imagetriu(medmat(diagidx,diagidx),50,180,flipud(jet),0);
    title(dataselected(vizidx(ii)));
    caxis([150,250]);
    cbar = colorbar;
    cbar.Label.String = 'nm';
    cbar.Label.FontSize = 15;
    set(gca,'XTick',xlabelspos,'XTickLabels',xlabels,'FontSize',10,...
        'XColor','k','YColor','k');
    xlabel('kb');
%     saveas(fh,[figDIR,'avg_medmat_',dataselected{vizidx(ii)},'.pdf']);
end

%% AVG -- subtraction of median distancemap 
% SizedFig(60,10);
fh = SizedFig(100,20);

% diagidx = [1,2,3,4,7,10,13,16,17,18,19];%1:19;
% diagidx = 1:13;
% diagidx = 4:16;
diagidx = 1:19;

vizidx = [1,2,3,4];
% vizidx = [7,8,10,11,13];
% vizidx = 1:length(medmats);%[1,3,4,5];
% vizidx = [2,3, 7,8];
% vizidx = [2,3, 7,8];
% vizidx = [9,10,11];
xlabels = {'–75','–30','–15','0','15','30','75'};
xlabelspos = [1,4,7,10,13,16,19];
for ii = 1:length(vizidx)
    subplot(1,length(vizidx),ii);
%     medmat = nanmean(medmats(vizidx(ii)).dmat - medmats(ctrlidx(vizidx(ii))).dmat,3);
    medmat = nanmean(medmats(vizidx(ii)).dmatfilt - medmats(ctrlidx(vizidx(ii))).dmatfilt,3);
    imagetriu(medmat(diagidx,diagidx),-50,50,flipud(bluewhiteredw0),0);
%     medmat = nanmean(medmats(vizidx(ii)).Rgfilt - medmats(ctrlidx(vizidx(ii))).Rgfilt,3);
%     imagetriu(medmat(diagidx,diagidx),-50,50,flipud(bluewhiteredw0),0);
    title(dataselected(vizidx(ii)));
    caxis([-50,50]);
    cbar = colorbar;
    cbar.Label.String = 'nm';
    cbar.Label.FontSize = 15;
    set(gca,'XTick',xlabelspos,'XTickLabels',xlabels,'FontSize',10);
    xlabel('kb');
%     saveas(fh,[figDIR,'subtraction_medmat_',dataselected{vizidx(ii)},'.pdf']);
end
%% local compaction by Rg Ratio; bar mean 

ii=4;
diagidx = 1:13;
medmatd = medmats(ii).Rg(diagidx,diagidx,:);
medmatc = medmats(ctrlidx(ii)).Rg(diagidx,diagidx,:);
localRgD = [];
localRgC = [];
for jj = 1:(size(medmatd,1)-4)
    localRgD(jj,:) = medmatd(jj,jj+4,:);
    localRgC(jj,:) = medmatc(jj,jj+4,:);
end
Rgratio = localRgD./localRgC;
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
ylim([0.86,1.12]);
% ylim([0.78,1.25]);
xlim([-30,30]);
xlabel('position [kb]');
ylabel('Rg ratio (treat/nodox)');
set(gca,'FontSize',15,'XColor','k','YColor','k');
box on;


%%
DIR = '/Users/tfuji/Documents/tempanalysis/_OnlyReporterPos/';
vizidx = [4];
ii = 1;
diagidx = 4:16;
medmat = nanmean(medmats(vizidx(ii)).dmat - medmats(ctrlidx(vizidx(ii))).dmat,3);
writematrix(medmat(diagidx,diagidx),[DIR,'dmatdiff_',dataselected{vizidx(ii)},'.csv']);
%%
diagidx = 4:16;
dmatintpall = [];
for ii = 1:length(medmats)
    tmpdmat = medmats(ii).dmat(diagidx,diagidx,:);
    dmatintpall = cat(3,dmatintpall,tmpdmat);
end
neighidx = triu(ones(size(dmatintpall,1),size(dmatintpall,2)),1)-triu(ones(size(dmatintpall,1),size(dmatintpall,2)),2);
neighdist = [];
for ii = 1:size(dmatintpall,3)
    tmp = dmatintpall(:,:,ii);
    neighdist = [neighdist;reshape(tmp(neighidx == 1),[],1)];
end
%% adjacent segment distance
SizedFig(20,20);
hold on;
xax = [-75,-60,-45,-30:5:30,45,60,75]';
xax = (xax(1:end-1)+xax(2:end))/2;
xlabels = {'–30','–20','–10','0','10','20','30'};

visidx = [2,3,4];%1:length(dataselected););
cmap = parula(length(visidx)+2);
cmap = cmap(2:end-1,:);
% cmap = [0.6 0.1 0.6; 0.7 0.4 0.1;0.1 0.6 0.7];
plot([-80,80],[0,0],'k--','LineWidth',1);    
for ii = 1:length(visidx)
    adjdist = nan(length(xax),5);
    medmat = medmats(visidx(ii)).dmat - medmats(ctrlidx(visidx(ii))).dmat;
    for jj = 1:(size(medmat,1)-1)
        adjdist(jj,:) = medmat(jj,jj+1,:);
    end
    avgdist = nanmean(adjdist,2);
    stddist = nanstd(adjdist,[],2) / sqrt(sum(~isnan(adjdist(1,:))));
    plot(xax,avgdist,'o-','Color',cmap(ii,:),'LineWidth',2);
    plot([xax,xax]',[avgdist-stddist,avgdist+stddist]','-','Color',cmap(ii,:),'LineWidth',1);
end
xlim([-80,80]);
xlim([-30,30]);
ylim([-20,20]);
set(gca, 'YDir','reverse','FontSize',10);
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
        
%% pick up upper triangle for dimentionality reduction
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
dfall = [];
for idxobj = 1:length(Objs)
    dfall = cat(2,dfall,Objs(idxobj).df);
end
idxrange = 4:16;
% idxrange = 8:14;
for ii = 1:length(dfall)
    allpairs = [];
    tmp = dfall(ii).dmatinterp(idxrange,idxrange,:);
    tmp_triu = ~triu(ones(size(tmp,1),size(tmp,2)))';
    for i = 1:size(tmp,3)
        tmpmat = tmp(:,:,i);
        allpairs = cat(1,allpairs,tmpmat(tmp_triu)');
    end
    dfall(ii).allpairs = allpairs;
    idx = 1:size(tmp,3);%randsample(size(tmp,3),267);
    dfall(ii).pairs = allpairs(idx,:);
end

pairs_allsample = [];
dmat_allsample = [];
coord_allsample = [];
datarange = 1:length(dfall);%[1,2,3,4,5,6,7,8,9];%

for ii = datarange%1:length(df)
    pairs_allsample = cat(1,pairs_allsample,dfall(ii).pairs);
    dmat_allsample = cat(3,dmat_allsample,dfall(ii).dmatinterp(idxrange,idxrange,:));
    coord_allsample = cat(3,coord_allsample,dfall(ii).coordinterp(idxrange,:,:));
end
%% display names of data
for ii = 1:length(dfall)
    disp([num2str(ii),':',dfall(ii).dname]);
end
dataselected = {'ctrl','difday1','difday2','difday3'};
ctrlidx = [1,1,1,1];
IdxByCond = repmat({[]},length(dataselected),1);
for jj = 1:length(dataselected)
    for ii = 1:length(dfall)
        if strcmp(dfall(ii).dname,dataselected{jj})
            IdxByCond{jj} = cat(2,IdxByCond{jj},ii);
        end
    end
end

%% k-means
rng('default');
kmnum = 16;
% kmnum = 6;
kmopt = kmeans(pairs_allsample,kmnum);

%
sizelist = [];
for ii = 1:length(datarange)
    sizelist = cat(1,sizelist,size(dfall(datarange(ii)).pairs,1));
end
idxst = cat(1,1,cumsum(sizelist(1:end-1))+1);
idxed = cumsum(sizelist);

for ii = 1:length(datarange)
    dfall(datarange(ii)).km = kmopt(idxst(ii):idxed(ii),:);
end

%% avg of each cluster
SizedFig(50,20);
for ii = 1:kmnum
    subplot(ceil(sqrt(kmnum)),ceil(sqrt(kmnum)),ii);
    tmp = dmat_allsample(:,:,kmopt == ii);
%     imagetriu(nanmedian(tmp,3),150,300,flipud(jet));
    imagetriu(nanmedian(tmp,3),100,400,flipud(jet));
    title(num2str(ii));
end

%% k_mean percentage by exp condition and bio rep
visidx = IdxByCond; % KRABnodox,irrversible,reactivated%
dnum = length(visidx);

kmean_tblall = nan(kmnum,dnum,5);
Rg_all = nan(1,dnum,5);
for ii = 1:dnum
tmpidx = visidx{ii};
tmp = [];
    for jj = 1:length(tmpidx)
        tbl = tabulate(dfall(tmpidx(jj)).km);
%         tmpRg = reshape(nanmedian(dfall(tmpidx(jj)).Rg(1,13,:)),[],1);
        idx = dfall(tmpidx(jj)).rep;
        kmean_tblall(:,ii,idx) = tbl(:,3);
%         Rg_all(:,ii,idx) = tmpRg;
    end
end
%% k_mean percentage bar graph / scatter
visidx = [1,2,3,4];
conditionnum = length(visidx);
cmap = parula(size(kmean_tblall,1)+2);
cmap = cmap(2:end-1,:);
wd = 0.2; 
ylims = [30,50];
SizedFig(50,30); 
for idxkm = 1:size(kmean_tblall,1)
    subplot(ceil(sqrt(size(kmean_tblall,1))),...
        ceil(sqrt(size(kmean_tblall,1))),idxkm);
    % idxkm = 4; ylims = [0,55];
    % idxkm = 1; ylims = [0,22];
    for i = 1:conditionnum
        kk = visidx(i);
        hold on;
        tmp = reshape(kmean_tblall(idxkm,kk,:),[],1);
        tmpref = reshape(kmean_tblall(idxkm,1,:),[],1);
        tmpref = tmpref(~isnan(tmp));
        tmp = tmp(~isnan(tmp));
        bar(i,nanmean(tmp),'FaceColor',cmap(idxkm,:),'LineWidth',2);
        xax = linspace(i-wd,i+wd,length(tmp));
        plot(xax,tmp,'ko','MarkerFaceColor','w','MarkerSize',8,...
            'LineWidth',2);
        [~,pv] = ttest(tmp,tmpref);
        disp(pv);
        box on;
    end
%     ylim(ylims);
    set(gca,'FontSize',20,'XTick',[],'XColor','k','YColor','k');
    % xlabel('after dox removal [days]');
    ylabel('memory fraction [%]');
    % legend(string(rec_days),'Location','eastoutside');
end

%% k_mean percentage bar graph / scatter
tmpdata = kmean_tblall(1,:,:);%+kmean_tblall(2,:,:);
visidx = [1,2,3,4];
conditionnum = length(visidx);
cmap = parula(size(tmpdata,1)+2);
cmap = cmap(2:end-1,:);
wd = 0.2; 
ylims = [30,50];
SizedFig(20,30); 
for i = 1:conditionnum
    kk = visidx(i);
    hold on;
    tmp = reshape(tmpdata(1,kk,:),[],1);
    tmpref = reshape(tmpdata(1,1,:),[],1);
    tmpref = tmpref(~isnan(tmp));
    tmp = tmp(~isnan(tmp));
    bar(i,nanmean(tmp),'FaceColor',cmap(1,:),'LineWidth',2);
    xax = linspace(i-wd,i+wd,length(tmp));
    plot(xax,tmp,'ko','MarkerFaceColor','w','MarkerSize',8,...
        'LineWidth',2);
    [~,pv] = ttest(tmp,tmpref);
    disp(pv);
    box on;
end
%     ylim(ylims);
set(gca,'FontSize',20,'XTick',[],'XColor','k','YColor','k');
% xlabel('after dox removal [days]');
ylabel('memory fraction [%]');
% legend(string(rec_days),'Location','eastoutside');

%% tSNE
rng('default');
tSNEout = tsne(pairs_allsample,'Perplexity',35,'Algorithm','barneshut','Distance','euclidean');

%
sizelist = [];
for ii = 1:length(datarange)
    sizelist = cat(1,sizelist,size(dfall(datarange(ii)).pairs,1));
end
idxst = cat(1,1,cumsum(sizelist(1:end-1))+1);
idxed = cumsum(sizelist);

for ii = 1:length(datarange)
    dfall(datarange(ii)).tSNEout = tSNEout(idxst(ii):idxed(ii),:);
end

%% plot t-sne for individual condition
% SizedFig(40,50);
SizedFig(30,10);
% visidx = [19,23,26,29,30,34];
% visidx = 1:9;
% visidx = 10:18;
visidx = [1,2,3,4];
% visidx = [10,11,12,33,34];
dnum = length(visidx);
for ii = 1:dnum
% subplot(ceil(sqrt(dnum)),ceil(sqrt(dnum)),ii);
subplot(1,dnum,ii);
dscatter(dfall(visidx(ii)).tSNEout(:,1),dfall(visidx(ii)).tSNEout(:,2));
xlabel('t-SNE 1');
ylabel('t-SNE 2');
title(dfall(visidx(ii)).dname);
axis image;
% axis([-50,50,-50,50]);
% axis([-60,60,-60,60]);
caxis([0.2*10^-4,0.5*10^-4]);
colormap(jet);

hold on;
box on;
end

%% plot t-sne for all condition; colored by k-mean
SizedFig(10,15);
tmptsne = [];
tmpkm = [];
for ii = 1:length(dfall)
    tmptsne = cat(1,tmptsne,dfall(ii).tSNEout);
    tmpkm   = cat(1,  tmpkm,dfall(ii).km);
end
cmap = parula(max(unique(tmpkm))+2);
cmap = cmap(2:end-1,:);
scatter(tmptsne(:,1),tmptsne(:,2),8,tmpkm,'filled');
colormap(gca,cmap);
xlabel('t-SNE 1');
ylabel('t-SNE 2');
axis image;
axis([-50,50,-50,50]);
set(gca,'XColor','k','YColor','k');

box on;
%% kernel density estimation
cmap = parula(256);
axrng = [-50,50];
% axrng = [-60,60];
figure;
for idx = 1:length(dfall)
    [bandwidth,densityD,X,Y]=kde2d(dfall(idx).tSNEout,32,[axrng(1),axrng(1)],[axrng(2),axrng(2)]);
    densityD = densityD/sum(densityD(:));
    % plot the data and the density estimate
    imagesc(densityD)
    colormap hot, hold on, alpha(.8)
    caxis([0,0.1*10^-3]);
    set(gca, 'color', 'blue');
    dfall(idx).kdemat = densityD;
end
%% average density tsne
SizedFig(100,20);
visidx = IdxByCond;%([1,2,3,4,5]); % KRABnodox,irrversible,reactivated%
dnum = length(visidx);

kdematall = [];
for ii = 1:dnum
subplot(1,dnum,ii);
tmpidx = visidx{ii};
tmp = [];
    for jj = 1:length(tmpidx)
%         tmp = cat(3,tmp,dfall(tmpidx(jj)).densemat);
        tmp = cat(3,tmp,dfall(tmpidx(jj)).kdemat);
        disp(dfall(tmpidx(jj)).dname);
    end
    imagesc(mean(tmp,3));
    kdematall = cat(3,kdematall,mean(tmp,3));
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    title(dfall(tmpidx(1)).dname);
    axis image; axis xy;
%     caxis([0.0*10^-5,0.8*10^-2]);
    caxis([1.0*10^-5,2.2*10^-3]);
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    hold on;
end
%%
corrmat = zeros(32,32);
for ii = 1:32
    for jj = 1:32
        tmp = corrcoef([2,1,3,4],reshape(kdematall(ii,jj,:),[],1));
        corrmat(ii,jj) = tmp(1,2);
    end
end
figure; imagesc(corrmat); axis image; axis xy;
%% average density tsne; subtraction

SizedFig(30,10);
SizedFig(100,20);
% visidx = ([1,2,3,4,5]); % KRABnodox,irrversible,reactivated % IdxByCond;%
visidx = 1:4; % KRABnodox,irrversible,reactivated % IdxByCond;%
dnum = length(visidx);

for ii = 1:dnum%size(kdematall,3)
    subplot(1,dnum,ii);
    tmp = kdematall(:,:,ii) - kdematall(:,:,ctrlidx(ii));
    imagesc(tmp);
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    title(dataselected(ii));
    axis image; axis xy;
    colormap(bluewhiteredw0);
    caxis([-1,1]*0.3*10^-2);
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    hold on;
    sum(sum(kdematall(:,:,ii).*log(kdematall(:,:,ii)./kdematall(:,:,ctrlidx(ii)))))
end

%% histogram of Rg
idx1 = 1;
idx2 = 4;
seg1 = 7; seg2 = 12;
SizedFig(20,20);
tmp1 = reshape(dfall(idx1).Rg(seg1,seg2,:),[],1);
histogram(tmp1,50,'BinLimits',[0,400],'normalization','probability');
hold on;
tmp2 = reshape(dfall(idx2).Rg(seg1,seg2,:),[],1);
histogram(tmp2,50,'BinLimits',[0,400],'normalization','probability');
ranksum(tmp1,tmp2)

%% histogram of Rg
idx1 = 1;
idx2 = 4;
seg1 = 7; seg2 = 12;
SizedFig(20,20);
tmp1 = log(reshape(dfall(idx1).Rg(seg1,seg2,:),[],1));
histogram(tmp1,60,'BinLimits',[3,6],'normalization','probability');
hold on;
tmp2 = log(reshape(dfall(idx2).Rg(seg1,seg2,:),[],1));
histogram(tmp2,60,'BinLimits',[3,6],'normalization','probability');
