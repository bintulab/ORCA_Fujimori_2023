%% directory for figure
figDIR = 'Z:\Taihei\project_ORCA\singleclones_ORCA\Analysis\fig1\';
%% re-aggregate data

dataselected = {'KRAB-nodox','KRAB-1day','KRAB-5days','reactivated','irreversible',...
    'KRAB147-nodox','KRAB147-1day','KRAB147-5days','KRAB150-nodox','KRAB150-1day','KRAB150-5days',...
    'HDAC4-nodox','HDAC4-1day','VP64-nodox','VP64-1day'};
ctrlidx = [1,1,1,1,1,6,6,6,9,9,9,12,12,14,14];

dfig1 = [];
dfig1.dmat = [];
dfig1.dmatNoIntp = [];
dfig1.dmatfilt = [];
dfig1.coordfilt = [];
dfig1.coordNoIntp = [];
dfig1.Rg = [];
dfig1.Rgfilt = [];
dfig1.rep = [];
dfig1.rep_dmatfilt = [];
dfig1.ctrlidx = [];
dfig1.ctrlidx_dmatfilt = [];
dfig1.dname = {};
dfig1.dname_dmatfilt = {};
for kk = 1:length(Objs)
   df1 = Objs(kk).df;
    drange = 1:length(df1);
    for ii = 1:length(drange)
        if sum(strcmp(df1(ii).dname,dataselected))
            dnum = size(df1(ii).dmatinterp,3);
            dfig1.dmat = cat(3,dfig1.dmat,df1(ii).dmatinterp);
            tmpRg = zeros(13,13,dnum);
            for jj = 4:15
                for mm = (jj+1):16
                    idxrange = jj:mm;
                    tmp = df1(ii).coordinterp(idxrange,:,:);
                    tmp = tmp - repmat(mean(tmp,1),[size(tmp,1), 1, 1]);
                    tmpRg(jj-3,mm-3,:) = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
                    tmpRg(mm-3,jj-3,:) = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
                end
            end
            dfig1.Rg = cat(3,dfig1.Rg,tmpRg);
            dfig1.dmatNoIntp = cat(3,dfig1.dmatNoIntp,df1(ii).dmatfilt50p);
            dfig1.coordNoIntp = cat(3,dfig1.coordNoIntp,df1(ii).coordfilt50p);
            dfig1.dmatfilt = cat(3,dfig1.dmatfilt,df1(ii).dmatfilt);
            dfig1.coordfilt = cat(3,dfig1.coordfilt,df1(ii).coordfilt);
            for jj = 1:dnum
                dfig1.rep = cat(1,dfig1.rep,kk);            
                dfig1.ctrlidx = cat(1,dfig1.ctrlidx,ctrlidx(strcmp(df1(ii).dname,dataselected)));            
                dfig1.dname = [dfig1.dname;df1(ii).dname];            
            end
            
            tmpRg = zeros(13,13,size(df1(ii).dmatfilt,3));
            for jj = 4:15
                for mm = (jj+1):16
                    idxrange = jj:mm;
                    tmp = df1(ii).coordfilt(idxrange,:,:);
                    tmp = tmp - repmat(mean(tmp,1),[size(tmp,1), 1, 1]);
                    tmpRg(jj-3,mm-3,:) = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
                    tmpRg(mm-3,jj-3,:) = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
                end
            end
            dfig1.Rgfilt = cat(3,dfig1.Rgfilt,tmpRg);
            for jj = 1:size(df1(ii).dmatfilt,3)
                dfig1.dname_dmatfilt = [dfig1.dname_dmatfilt;df1(ii).dname];            
                dfig1.rep_dmatfilt = cat(1,dfig1.rep_dmatfilt,kk);            
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
    pmat = [];
    ksmat = [];
    nposmat = [];
    for jj = 1:length(Objs)
        tmp = dfig1.dmat(:,:,strcmp(dfig1.dname,dataselected(ii)) & dfig1.rep == jj);
        tmp = nanmedian(tmp,3);
        medmat = cat(3,medmat,tmp);
        % example plot
%         subplot(1,length(Objs),jj);
%         if ~isnan(mean(tmp(:)))
%             imagetriu(tmp(diagidx,diagidx),200,300,flipud(jet),0);
%             title(dataselected(ii));
%             caxis([200,300]);
%             cbar = colorbar;
%             cbar.Label.String = 'nm';
%             cbar.Label.FontSize = 15;
%             set(gca,'XTick',1:2:13,'XTickLabels',xlabels,'FontSize',10);
%         end
        tmp = dfig1.dmatfilt(:,:,strcmp(dfig1.dname_dmatfilt,dataselected(ii)) & dfig1.rep_dmatfilt == jj);
        npostmp = sum(~isnan(tmp),3)/size(tmp,3);
        nposmat = cat(3,nposmat,npostmp);
        tmp = nanmedian(tmp,3);
        medmat2 = cat(3,medmat2,tmp);
        tmp = dfig1.Rg(:,:,strcmp(dfig1.dname,dataselected(ii)) & dfig1.rep == jj);
        tmp = nanmedian(tmp,3);
        medmat3 = cat(3,medmat3,tmp);
        tmp = dfig1.Rgfilt(:,:,strcmp(dfig1.dname_dmatfilt,dataselected(ii)) & dfig1.rep_dmatfilt == jj);
        tmp = nanmedian(tmp,3);
        medmat4 = cat(3,medmat4,tmp);

        % calc p-val
        tmp1 = dfig1.dmatfilt(:,:,strcmp(dfig1.dname_dmatfilt,dataselected(ii)) & dfig1.rep_dmatfilt == jj);
        tmp2 = dfig1.dmatfilt(:,:,strcmp(dfig1.dname_dmatfilt,dataselected(ctrlidx(ii))) & dfig1.rep_dmatfilt == jj);
        pmattmp  = nan(size(tmp1,1),size(tmp1,2));
        ksmattmp = nan(size(tmp1,1),size(tmp1,2));
        if size(tmp1,3) ~= 0
            for mm = 1:(size(tmp1,1)-1)
                for nn = (mm+1):size(tmp1,2)
                    x1 = reshape(tmp1(mm,nn,:),[],1);
                    x2 = reshape(tmp2(mm,nn,:),[],1);
                    rsum = ranksum(x1,x2);
                    [~,kstst] = kstest2(x1,x2);
                    pmattmp(mm,nn) = rsum;
                    pmattmp(nn,mm) = rsum;
                    ksmattmp(mm,nn) = kstst;
                    ksmattmp(nn,mm) = kstst;
                end
            end
        end
        pmat = cat(3,pmat,pmattmp);
        ksmat = cat(3,ksmat,ksmattmp);
    end
    medmats(ii).dmat = medmat;
    medmats(ii).dmatfilt = medmat2(1:end-1,1:end-1,:);
    medmats(ii).Rg = medmat3;
    medmats(ii).Rgfilt = medmat4;
    medmats(ii).pmatfilt = pmat;
    medmats(ii).ksmatfilt = ksmat;
    medmats(ii).nposmat = nposmat;
end
%% reproducibility check; rep1 med vs rep2 med
rep1 = [];
rep2 = [];
for ii = 1:length(medmats)
    tmp = medmats(ii).dmatfilt;
    nseg = size(tmp,1);
    triuidx = triu(ones(nseg)*1,1);
    for jj = 1:size(tmp,3)
        for kk = (jj+1):(size(tmp,3)-1)
            if ~isnan(tmp(1,1,jj)) && ~isnan(tmp(1,1,kk))
                tmp1 = tmp(:,:,jj);
                tmp2 = tmp(:,:,kk);
                rep1 = [rep1;tmp1(triuidx == 1)];
                rep2 = [rep2;tmp2(triuidx == 1)];
            end
        end
    end
end
%
SizedFig(15,20);
hold on;
% scatter(rep1,rep2,1,'b');
dscatter(rep1,rep2);
box on;
set(gca,'XColor','k','XColor','k');
xlim([100,500]);
ylim([100,500]);
refline(1,0);
corrcoef(rep1,rep2)
xlabel('rep1 median distance [nm]');
ylabel('rep2 median distance [nm]');


%% subtracted median distancemap for each biorep
% SizedFig(40,10);
diagidx = 4:16;
% diagidx = 1:13;
diagidx = 1:19;
xlabels = {'–75','–30','–15','0','15','30','75'};
xlabelspos = [1,4,7,10,13,16,19];
for ii = [5]%[2,3,13]%1:length(dataselected)%

    SizedFig(80,30);
    for jj = 7:length(Objs)
        subplot(1,length(Objs),jj);
%         tmp1 = dfig1.dmat(:,:,strcmp(dfig1.dname,dataselected(ii)) & dfig1.rep == jj);
%         tmp1 = nanmedian(tmp1,3);
%         tmp2 = dfig1.dmat(:,:,strcmp(dfig1.dname,dataselected(ctrlidx(ii))) & dfig1.rep == jj);
%         tmp2 = nanmedian(tmp2,3);
        tmp1 = dfig1.dmatfilt(:,:,strcmp(dfig1.dname_dmatfilt,dataselected(ii)) & dfig1.rep_dmatfilt == jj);
        tmp1 = nanmedian(tmp1,3);
        tmp2 = dfig1.dmatfilt(:,:,strcmp(dfig1.dname_dmatfilt,dataselected(ctrlidx(ii))) & dfig1.rep_dmatfilt == jj);
        tmp2 = nanmedian(tmp2,3);
        tmp = tmp1 - tmp2;
        if ~isnan(mean(tmp(:)))
            imagetriu(rot90(tmp(diagidx,diagidx),2),-100,100,flipud(bluewhiteredw0),0);
            caxis([-100,100]);
%             imagetriu(rot90(tmp1(diagidx,diagidx),2),250,350,flipud(jet),0);
%             caxis([250,350]);
            title(dataselected(ii));
            cbar = colorbar;
            cbar.Label.String = 'nm';
            cbar.Label.FontSize = 15;
            set(gca,'XTick',xlabelspos,'XTickLabels',xlabels,'FontSize',10,...
                'XColor','k','YColor','k');
        end
    end
end
%% AVG -- median distancemap 
% SizedFig(60,10);
fh = SizedFig(100,20);

diagidx = 4:16;
xlabels = {'–30','–15','0','15','30'};
xlabelspos = [1,4,7,10,13];

% diagidx = 1:19;
% xlabels = {'–75','–30','–15','0','15','30','75'};
% xlabelspos = [1,4,7,10,13,16,19];

vizidx = [1,3,4,5];%1:length(medmats);%
% vizidx = [9,10,11];%1:length(medmats);%
for ii = 1:length(vizidx)
    subplot(1,length(vizidx),ii);
    medmat = nanmean(medmats(vizidx(ii)).dmatfilt,3);
    imagetriu(rot90(medmat(diagidx,diagidx),2),230,350,flipud(jet),0);
%     medmat = nanmean(medmats(vizidx(ii)).Rg,3);
%     imagetriu(medmat(diagidx,diagidx),50,180,flipud(jet),0);
    title(dataselected(vizidx(ii)));
    caxis([230,350]);
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
% fh = SizedFig(100,20);

% diagidx = [1,2,3,4,7,10,13,16,17,18,19];%1:19;
% diagidx = 1:13;

diagidx = 4:16;
xlabels = {'–30','–15','0','15','30'};
xlabelspos = [1,4,7,10,13];

% diagidx = 1:19;
% xlabels = {'–75','–30','–15','0','15','30','75'};
% xlabelspos = [1,4,7,10,13,16,19];

vizidx = [1,3,4,5];%[3];%
% vizidx = [2,7,8,10,11,13];
% vizidx = 1:length(medmats);%[1,3,4,5];
% vizidx = [2,3, 7,8];
% vizidx = [2,3, 7,8]
% vizidx = [9,10,11];
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
%% histogram of subtracted median distance
cmap = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980];
ii=5;
medmat = nanmean(medmats(ii).dmatfilt - medmats(ctrlidx(ii)).dmatfilt,3);
meddif_wreporter = [];
meddif_woreporter = [];
meddif_all = [];
for ii = 1:(size(medmat,1)-1)
    for jj = (ii+1):size(medmat,1)
        if ii <= 10 && jj >= 10
            meddif_wreporter = [meddif_wreporter;medmat(ii,jj)];
        else
            meddif_woreporter = [meddif_woreporter;medmat(ii,jj)];
        end
        meddif_all = [meddif_all;medmat(ii,jj)];
    end
end
SizedFig(20,25); 
hold on;
h1 = histogram(meddif_woreporter,12,'BinLimits',[-70,30],'normalization','probability',...
    'EdgeAlpha',0.0,'FaceColor',cmap(1,:));
stairs([h1.BinEdges,h1.BinEdges(end)],[h1.Values,h1.Values(end),0],'k-');
h2 = histogram(meddif_wreporter,12,'BinLimits',[-70,30],'normalization','probability',...
    'EdgeAlpha',0.0,'FaceColor',cmap(2,:));
stairs([h2.BinEdges,h2.BinEdges(end)],[h2.Values,h2.Values(end),0],'k-');
plot([0,0],[0,0.35],'k--');

set(gca,'XColor','k','YColor','k');
box on;
xlabel('difference in median');
ylabel('probability');

%% local compaction by Rg; bar mean,  xreversed
cmap1 = [0, 0.4470, 0.7410];
cmap2 = [0.8500, 0.3250, 0.0980];
{'KRAB-nodox','KRAB-1day','KRAB-5days','reactivated','irreversible',...
    'KRAB147-nodox','KRAB147-1day','KRAB147-5days','KRAB150-nodox','KRAB150-1day','KRAB150-5days',...
    'HDAC4-nodox','HDAC4-1day','VP64-nodox','VP64-1day'};
ii=3;
diagidx = 1:13;
medmatd = medmats(ii).Rg(diagidx,diagidx,:);
medmatc = medmats(ctrlidx(ii)).Rg(diagidx,diagidx,:);
localRgD = [];
localRgC = [];
for jj = 1:(size(medmatd,1)-4)
    localRgD(jj,:) = medmatd(jj,jj+4,:);
    localRgC(jj,:) = medmatc(jj,jj+4,:);
end
nonnanidx = ~isnan(localRgD(1,:));
localRgD = flipud(localRgD);
localRgD = localRgD(:,nonnanidx);
avgRgRatio = nanmean(localRgD,2);
stdRgRatio = nanstd(localRgD,[],2) / sqrt(sum(~isnan(localRgD(1,:))));
localRgC = flipud(localRgC);
localRgC = localRgC(:,nonnanidx);
avglocalRgC = nanmean(localRgC,2);
stddistlocalRgC = nanstd(localRgC,[],2) / sqrt(sum(~isnan(localRgC(1,:))));

SizedFig(30,25);
hold on;
xax = linspace(-floor(size(localRgD,1)/2)*5,floor(size(localRgD,1)/2)*5,size(localRgD,1));
% plot(repmat(xax,[2,1]),[avglocalRgC-stddistlocalRgC,avglocalRgC+stddistlocalRgC]','k-','linewidth',1);
pgon = polyshape([xax,fliplr(xax)],[avglocalRgC-stddistlocalRgC;flipud(avglocalRgC+stddistlocalRgC)]');
plot(pgon,'EdgeColor','none','FaceColor',cmap1,'FaceAlpha',0.5);
plot(xax,avglocalRgC,'.-','linewidth',2,'Color',cmap1,'MarkerSize',20);

pgon = polyshape([xax,fliplr(xax)],[avgRgRatio-stdRgRatio;flipud(avgRgRatio+stdRgRatio)]');
plot(pgon,'EdgeColor','none','FaceColor',cmap2,'FaceAlpha',0.5);
plot(xax,avgRgRatio,'.-','linewidth',2,'Color',cmap2,'MarkerSize',20);

xlim([-20,20]);
ylim([118,170]);
xlabel('position [kb]');
ylabel('Local radius of gyration (±10kb) [nm]');
set(gca,'FontSize',15,'XColor','k','YColor','k');
box on;
%%
disp('---');
for ii = 1:size(localRgC,1)
    [~,pval] = ttest(localRgC(ii,:),localRgD(ii,:));
    disp(pval);
end

%% local compaction by Rg Ratio; bar mean, xreversed
{'KRAB-nodox','KRAB-1day','KRAB-5days','reactivated','irreversible',...
    'KRAB147-nodox','KRAB147-1day','KRAB147-5days','KRAB150-nodox','KRAB150-1day','KRAB150-5days',...
    'HDAC4-nodox','HDAC4-1day','VP64-nodox','VP64-1day'};
ii=5;
diagidx = 1:13;
medmatd = medmats(ii).Rg(diagidx,diagidx,:);
medmatc = medmats(ctrlidx(ii)).Rg(diagidx,diagidx,:);
localRgD = [];
localRgC = [];
for jj = 1:(size(medmatd,1)-4)
    localRgD(jj,:) = medmatd(jj,jj+4,:);
    localRgC(jj,:) = medmatc(jj,jj+4,:);
end
%
% localRgD = localRgD(:,[4,5,6]);
% localRgC = localRgC(:,[4,5,6]);
nonnanidx = ~isnan(localRgD(1,:));
localRgD = localRgD(:,nonnanidx);
localRgC = localRgC(:,nonnanidx);
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
ylim([0.87,1.12]);
% ylim([0.78,1.25]);
xlim([-30,30]);
xlabel('position [kb]');
ylabel('Rg ratio (treat/nodox)');
set(gca,'FontSize',15,'XColor','k','YColor','k');
box on;

%% AVG -- rank sum pval map xreversed
% SizedFig(60,10);
fh = SizedFig(100,20);

% diagidx = [1,2,3,4,7,10,13,16,17,18,19];%1:19;
% diagidx = 1:13;
% diagidx = 4:16;
diagidx = 1:19;

vizidx = [1,3,4,5];
% vizidx = [7,8,10,11,13];
% vizidx = 1:length(medmats);%[1,3,4,5];
% vizidx = [2,3, 7,8];
% vizidx = [2,3, 7,8];
% vizidx = [9,10,11];
xlabels = {'–75','–30','–15','0','15','30','75'};
xlabelspos = [1,4,7,10,13,16,19];
for ii = 1:length(vizidx)
    subplot(1,length(vizidx),ii);
    tmp = medmats(vizidx(ii)).pmatfilt;
    for jj = 1:size(tmp,1)
        tmp(jj,jj) = 1;
    end
    medmat = nanmean(log10(tmp),3);
%     medmat = nanmin(medmats(vizidx(ii)).pmatfilt,[],3);
    imagetriu(rot90(medmat(diagidx,diagidx),2),-1.2,-0.2,flipud(hot),0);
    title(dataselected(vizidx(ii)));
    caxis([-1.2,-0.2]);
    cbar = colorbar;
    cbar.Label.String = 'log10 p';
    cbar.Label.FontSize = 15;
    set(gca,'XTick',xlabelspos,'XTickLabels',xlabels,'FontSize',10);
    xlabel('kb');
%     saveas(fh,[figDIR,'subtraction_medmat_',dataselected{vizidx(ii)},'.pdf']);
end



%% representative 3D structures; closest avg median map
dataselected = {'KRAB-nodox','KRAB-1day','KRAB-5days','reactivated','irreversible',...
    'KRAB147-nodox','KRAB147-1day','KRAB147-5days','KRAB150-nodox','KRAB150-1day','KRAB150-5days',...
    'HDAC4-nodox','HDAC4-1day','VP64-nodox','VP64-1day'};
dataidx = 1;
azseries = [53,53,53+90];%+180;
dataidx = 3;
azseries = [53,53,53+90];%+180;
dataidx = 4;
azseries = [53,53-90,53-90];%+180;
dataidx = 5;
azseries = [53,53,53];%+180;
% tmpcoord = dfig1.coordfilt(1:(end-1),:,strcmp(dfig1.dname_dmatfilt,dataselected(dataidx)));
% endpos = ~isnan(tmpcoord(1,1,:)) & ~isnan(tmpcoord(19,1,:));
% tmpdmat  = dfig1.dmatfilt(1:(end-1),1:(end-1),strcmp(dfig1.dname_dmatfilt,dataselected(dataidx)));
% tmpcoord = tmpcoord(:,:,endpos);
% tmpdmat = tmpdmat(:,:,endpos);
drange = 4:16; % 1:19;
tmpcoord = dfig1.coordNoIntp(drange,:,strcmp(dfig1.dname,dataselected(dataidx)));
endpos = ~isnan(tmpcoord(1,1,:)) & ~isnan(tmpcoord(end,1,:));
% tmpdmat  = dfig1.dmatNoIntp(drange,drange,strcmp(dfig1.dname,dataselected(dataidx)));
tmpdmat  = dfig1.dmat(drange,drange,strcmp(dfig1.dname,dataselected(dataidx)));
tmpcoord = tmpcoord(:,:,endpos);
tmpdmat = tmpdmat(:,:,endpos);
%
refmat = nanmean(medmats(dataidx).dmatfilt(drange,drange,:),3);
difffrommed = reshape(nansum(nansum((tmpdmat - refmat).^2,1),2),[],1);
% difffrommed = reshape(sum(sum(abs(dmat1 - refmat),1),2),[],1);
idxclosest = find(difffrommed == min(difffrommed))
[~,minidx] = sort(difffrommed,'ascend');

elseries = [1,1,1]*20;
for jj = 1:3
    fh = SizedFig(20,30); 
%     subplot(1,3,jj);
    minidx(jj)
    repcoord = tmpcoord(:,:,minidx(jj));

    % visualize the 3D structure
    cntr = (drange(end) - drange(1))/2 + 1;
    az = azseries(jj);    el = elseries(jj);

    tmp = repcoord; 
    xangle = 0; yangle = 0; zangle = 0;
    % tmp = repcoordkrab; 
    % xangle = 0; yangle = 0; zangle = 0;
    cmapust = [255,0,0; 255,0,0;
        255,131,0; 255,131,0;
        248,255,0; 248,255,0;
    117,255,0; 117,255,0;
    73,255,0; 73,255,0;
    29,255,0; 29,255,0;
    0,255,15; 0,255,15;
    0,255,58; 0,255,58;
    0,255,102; 0,255,102]/255;
    ustidx = zeros(size(cmapust,1),1);
    ustidx(1:2:(cntr-1)*2) = flipud(~isnan(tmp(1:(cntr-1),1)));
    ustidx(2:2:(cntr-1)*2) = flipud(~isnan(tmp(1:(cntr-1),1)));
    ustidx = flipud(ustidx);
    cmaprep = [220, 165, 104]/255;
    cmapdst = [0,255,233; 0,255,233;
    0,233,255; 0,233,255;
    0,189,255; 0,189,255;
    0,146,255; 0,146,255;
    0,102,255; 0,102,255;
    0,58,255; 0,58,255;
        73,0,255; 73,0,255;
        204,0,255; 204,0,255;
        255,0,175; 255,0,175]/255;
    dstidx = zeros(size(cmapust,1),1);
    dstidx(1:2:(cntr-1)*2) = ~isnan(tmp((cntr+1):end,1));
    dstidx(2:2:(cntr-1)*2) = ~isnan(tmp((cntr+1):end,1));
    cmap = cat(1,cmapust(ustidx == 1,:),cmaprep,cmapdst(dstidx == 1,:));
    boxhw = 330;    

    tmp = tmp(~isnan(tmp(:,1)),:); 
    polycenter = (max(tmp,[],1) + min(tmp,[],1))/2;% mean(tmp,1);
    tmp = tmp - polycenter;
    % tmp = (rotz(zangle)*roty(yangle)*rotx(xangle)*tmp')';
%     tmp(:,3) = tmp(:,3)+10; % to frame the figure in 400x400x400 box. just for aesthetics
    title(dataselected{dataidx});    
    minmax = [min(tmp,[],1) - 0.3*abs(max(tmp,[],1) - min(tmp,[],1));...
        max(tmp,[],1) + 0.3*abs(max(tmp,[],1) - min(tmp,[],1))];

    if sum(abs(minmax(1,:) - minmax(2,:))) ~= 0 && ~isnan(sum(abs(minmax(1,:) - minmax(2,:))))
        tmpinterpx = [];
        tmpinterpy = [];
        tmpinterpz = [];
        tmpinterpx = cat(1,tmpinterpx,...
            interp1(1:size(tmp,1),tmp(:,1),1:0.1:size(tmp,1),'spline')');
        tmpinterpy = cat(1,tmpinterpy,...
            interp1(1:size(tmp,1),tmp(:,2),1:0.1:size(tmp,1),'spline')');
        tmpinterpz = cat(1,tmpinterpz,...
            interp1(1:size(tmp,1),tmp(:,3),1:0.1:size(tmp,1),'spline')');
        xbound = [-boxhw,boxhw];
        ybound = [-boxhw,boxhw];
        zbound = [-boxhw,boxhw];
        radius = 30;
        thickness = 10;
%         clf;%
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
                    thickness,cmap(i*2,:));
                % --- optimize axis --- %
                set(h, 'FaceLighting','phong','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
                set(h,'EdgeAlpha',0);
                material shiny; 
            end
        end

        % --- optimize axis --- %
        axis equal;
        grid off;
        set(gca,'Projection','perspective','Box','off','BoxStyle','full',...
            'FontSize',10)
        xlim(xbound); ylim(ybound); zlim(zbound);
%         xlabel('X [nm]');
%         ylabel('Y [nm]');
%         zlabel('Z [nm]');
        set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],...
            'XColor','k','YColor','k','XTick',[]);
    %         view(az(21),el);
        view(az,el);
    %         view([7.9927 28.5909]);
        camlight('right');
        cntr = 0;
        set(gca,'ZColor','none','XColor','none','YColor','none',...
            'LineWidth',2,'ZTick',[-200,0,200]);%,...
    %         'YGrid','off','XGrid','off','ZGrid','off');
    end
    figDIR = '/Users/tfuji/Documents/tempanalysis/_figures/fig1/rep_structure/';
    saveas(fh,[figDIR,'rep_3d_',dataselected{dataidx},'_rep_',num2str(jj),'.pdf']);
    
end
%%

rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)];
roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)];
rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1]; 


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