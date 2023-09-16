
%% directory for figure
figDIR = 'Z:\Taihei\project_ORCA\singleclones_ORCA\Analysis\fig1\';
%% re-aggregate data

dataselected = {'KRAB-nodox','KRAB-5days','reactivated','irreversible'};
ctrlidx = [1,1,1,1];

dataselected = {'KRAB-nodox','KRAB-1day','KRAB-5days','reactivated','irreversible',...
    'KRAB147-nodox','KRAB147-1day','KRAB147-5days','KRAB150-nodox','KRAB150-1day','KRAB150-5days',...
    'HDAC4-nodox','HDAC4-1day'};
ctrlidx = [1,1,1,1,1,6,6,6,9,9,9,12,12];

dfig4 = [];
dfig4.dmat = [];
dfig4.dmatNoIntp = [];
dfig4.dmatfilt = [];
dfig4.Rg = [];
dfig4.rep = [];
dfig4.rep_dmatfilt = [];
dfig4.ctrlidx = [];
dfig4.ctrlidx_dmatfilt = [];
dfig4.dname = {};
dfig4.dname_dmatfilt = {};
for kk = 1:2%length(Objs)
   df1 = Objs(kk).df;
    drange = 1:length(df1);
    for ii = 1:length(drange)
        if sum(strcmp(df1(ii).dname,dataselected))
            dnum = size(df1(ii).dmatinterp,3);
            dfig4.dmat = cat(3,dfig4.dmat,df1(ii).dmatinterp);
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
            dfig4.Rg = cat(3,dfig4.Rg,tmpRg);
            dfig4.dmatNoIntp = cat(3,dfig4.dmatNoIntp,df1(ii).dmatfilt50p);
            dfig4.dmatfilt = cat(3,dfig4.dmatfilt,df1(ii).dmatfilt);
            for jj = 1:dnum
                dfig4.rep = cat(1,dfig4.rep,kk);            
                dfig4.ctrlidx = cat(1,dfig4.ctrlidx,ctrlidx(strcmp(df1(ii).dname,dataselected)));            
                dfig4.dname = [dfig4.dname;df1(ii).dname];            
            end
            for jj = 1:size(df1(ii).dmatfilt,3)
                dfig4.dname_dmatfilt = [dfig4.dname_dmatfilt;df1(ii).dname];            
                dfig4.rep_dmatfilt = cat(1,dfig4.rep_dmatfilt,kk);            
            end
        end
    end
end

%% -----------------------------------------------------------------
%% --------- calculate median distance map for each bioreps --------
%% -----------------------------------------------------------------
%% median distancemap 
% SizedFig(40,10);
diagidx = 4:16;
% diagidx = 1:13;
% diagidx = 1:19;
xlabels = {'–30','–20','–10','0','10','20','30'};
medmats = [];
medmats.dmat = [];
medmats.Rg = [];
medmats.dmatNoIntp = [];
for ii = 1:length(dataselected)
%     SizedFig(40,10);
%     subplot(1,length(dataselected),ii);
    medmat = [];
    medmat2 = [];
    medmat3 = [];
    medmat4 = [];
    for jj = 1:2
        tmp = dfig4.dmat(:,:,strcmp(dfig4.dname,dataselected(ii)) & dfig4.rep == jj);
        tmp = nanmedian(tmp,3);
        medmat = cat(3,medmat,tmp);
%         subplot(1,2,jj);
%         if ~isnan(mean(tmp(:)))
%             imagetriu(tmp(diagidx,diagidx),200,300,flipud(jet),0);
%             title(dataselected(ii));
%             caxis([200,300]);
%             cbar = colorbar;
%             cbar.Label.String = 'nm';
%             cbar.Label.FontSize = 15;
%             set(gca,'XTick',1:2:13,'XTickLabels',xlabels,'FontSize',10);
%         end
        tmp = dfig4.dmatfilt(:,:,strcmp(dfig4.dname_dmatfilt,dataselected(ii)) & dfig4.rep_dmatfilt == jj);
        tmp = nanmedian(tmp,3);
        medmat2 = cat(3,medmat2,tmp);
        tmp = dfig4.Rg(:,:,strcmp(dfig4.dname,dataselected(ii)) & dfig4.rep == jj);
        tmp = nanmedian(tmp,3);
        medmat3 = cat(3,medmat3,tmp);

        % pval
        tmp1 = dfig4.Rg(1,end,strcmp(dfig4.dname,dataselected(ii)) & dfig4.rep == jj);
        tmp2 = dfig4.Rg(1,end,strcmp(dfig4.dname,dataselected(1)) & dfig4.rep == jj);
        if ~isempty(tmp1)
            tmp = ranksum(reshape(tmp1,[],1),reshape(tmp2,[],1));
            medmat4 = cat(1,medmat4,tmp);
        end
    end
    medmats(ii).dmat = medmat;
    medmats(ii).dmatfilt = medmat2(1:end-1,1:end-1,:);
    medmats(ii).Rg = medmat3;
    medmats(ii).Rg_p_vsKnodox = medmat4;
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
        tmp = dfig4.dmat(:,:,strcmp(dfig4.dname,dataselected(ii)) & dfig4.rep == jj);
        tmp = nanmedian(tmp,3);
        tmp2 = dfig4.dmat(:,:,strcmp(dfig4.dname,dataselected(ctrlidx(ii))) & dfig4.rep == jj);
        tmp2 = nanmedian(tmp2,3);
%         tmp = dfig4.dmatfilt(:,:,strcmp(dfig4.dname_dmatfilt,dataselected(ii)) & dfig4.rep_dmatfilt == jj);
%         tmp = nanmedian(tmp,3);
%         tmp2 = dfig4.dmatfilt(:,:,strcmp(dfig4.dname_dmatfilt,dataselected(ctrlidx(ii))) & dfig4.rep_dmatfilt == jj);
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

vizidx = [1,2,3,4,5];%1:length(medmats);%
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

vizidx = [1,3,4,5];
% vizidx = [7,8,10,11,13];
% vizidx = 1:length(medmats);%[1,3,4,5];
% vizidx = [2,3, 7,8];
vizidx = [2,3, 7,8];
diagidx = 1:19;
% xlabels = {'–75','–30','–15','0','15','30','75'};
% xlabelspos = [1,4,7,10,13,16,19];
% 
diagidx = 4:16;
xlabels = {'–30','–15','0','15','30'};
xlabelspos = [1,4,7,10,13];
for ii = 1:length(vizidx)
    subplot(1,length(vizidx),ii);
    medmat = nanmean(medmats(vizidx(ii)).dmatfilt - medmats(ctrlidx(vizidx(ii))).dmatfilt,3);
%     medmat = nanmean(medmats(vizidx(ii)).dmatfilt - medmats(1).dmatfilt,3);
    imagetriu(rot90(medmat(diagidx,diagidx),2),-100,100,flipud(bluewhiteredw0),0);
%     medmat = nanmean(medmats(vizidx(ii)).Rgfilt - medmats(ctrlidx(vizidx(ii))).Rgfilt,3);
%     imagetriu(medmat(diagidx,diagidx),-50,50,flipud(bluewhiteredw0),0);
    title(dataselected(vizidx(ii)));
    caxis([-100,100]);
    cbar = colorbar;
    cbar.Label.String = 'nm';
    cbar.Label.FontSize = 15;
    set(gca,'XTick',xlabelspos,'XTickLabels',xlabels,'FontSize',10);
    set(gca,'XColor','k','YColor','k');
    xlabel('kb');
%     saveas(fh,[figDIR,'subtraction_medmat_',dataselected{vizidx(ii)},'.pdf']);
end
%% local compaction by Rg Ratio; bar mean, xreversed
{'KRAB-nodox','KRAB-1day','KRAB-5days','reactivated','irreversible',...
    'KRAB147-nodox','KRAB147-1day','KRAB147-5days','KRAB150-nodox','KRAB150-1day','KRAB150-5days',...
    'HDAC4-nodox','HDAC4-1day','VP64-nodox','VP64-1day'};
ii=8;
diagidx = 1:13;
medmatd = medmats(ii).Rg(diagidx,diagidx,:);
% medmatc = medmats(ctrlidx(ii)).Rg(diagidx,diagidx,:);
medmatc = medmats(1).Rg(diagidx,diagidx,:);
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
ylim([0.86,1.12]);
ylim([0.78,1.25]);
xlim([-30,30]);
xlabel('position [kb]');
ylabel('Rg ratio (treat/nodox)');
set(gca,'FontSize',15,'XColor','k','YColor','k');
box on;

%% median distancemap 
% SizedFig(40,10);
diagidx = 4:16;
% diagidx = 1:13;
% diagidx = 1:19;
xlabels = {'–30','–20','–10','0','10','20','30'};
medmats = [];
medmats.dmat = [];
medmats.Rg = [];
medmats.dmatNoIntp = [];
for ii = 1:length(dataselected)
    SizedFig(40,10);
    subplot(1,length(dataselected),ii);
    medmat = [];
    medmat2 = [];
    medmat3 = [];
    for jj = 1:5
        subplot(1,5,jj);
        tmp = dfig4.dmat(:,:,strcmp(dfig4.dname,dataselected(ii)) & dfig4.rep == jj);
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
        tmp = dfig4.dmatfilt(:,:,strcmp(dfig4.dname_dmatfilt,dataselected(ii)) & dfig4.rep_dmatfilt == jj);
        tmp = nanmedian(tmp,3);
        medmat2 = cat(3,medmat2,tmp);
        tmp = dfig4.Rg(:,:,strcmp(dfig4.dname,dataselected(ii)) & dfig4.rep == jj);
        tmp = nanmedian(tmp,3);
        medmat3 = cat(3,medmat3,tmp);
    end
    medmats(ii).dmat = medmat;
    medmats(ii).dmatfilt = medmat2(1:end-1,1:end-1,:);
    medmats(ii).Rg = medmat3;
end


%% rank sum
visidx = [6,9,10,11,7,8,2,3];
SizedFig(30,30); 
for i = 1:length(visidx)
    hold on;
    bar(i,mean(log10(medmats(visidx(i)).Rg_p_vsKnodox)),'FaceColor','b','LineWidth',2);
    plot([i-0.2,i+0.2],log10(medmats(visidx(i)).Rg_p_vsKnodox),'ko','MarkerFaceColor','w','MarkerSize',8,...
        'LineWidth',2);
end
% ylim([0,100]);
set(gca,'FontSize',10,'XTick',1:length(visidx),'XTickLabel',dataselected(visidx),...
    'XColor','k','YColor','k');
ylabel('log10 p');
xtickangle(-45);
% legend(string(rec_days),'Location','eastoutside');
% xlabel('after dox removal [days]');
box on;