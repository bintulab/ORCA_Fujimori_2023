
%% directory for figure
figDIR = '/Users/tfuji/Documents/tempanalysis/_figures_v3/fig2/';
%% re-aggregate data
dataselected = {'KRAB-nodox','KRAB-1day','KRAB-5days','irreversible','reactivated',...
    'HDAC4-nodox','HDAC4-1day','HDAC4-5days','VP64-nodox','VP64-1day',...
    'DNMT3B-nodox','DNMT3B-irreversible',...
    'KRAB147-nodox','KRAB147-1day','KRAB147-5days',...
    'KRAB150-nodox','KRAB150-1day','KRAB150-5days'};
ctrlidx = [1,1,1,1,1,6,6,6,9,9,11,11,13,13,13,16,16,16];

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
    for jj = 1:length(Objs)
%         subplot(1,length(Objs),jj);
        tmp = dfig1.dmat(:,:,strcmp(dfig1.dname,dataselected(ii)) & dfig1.rep == jj);
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
        tmp = dfig1.dmatfilt(:,:,strcmp(dfig1.dname_dmatfilt,dataselected(ii)) & dfig1.rep_dmatfilt == jj);
        tmp = nanmedian(tmp,3);
        medmat2 = cat(3,medmat2,tmp);
        tmp = dfig1.Rg(:,:,strcmp(dfig1.dname,dataselected(ii)) & dfig1.rep == jj);
        tmp = nanmedian(tmp,3);
        medmat3 = cat(3,medmat3,tmp);
        tmp = dfig1.Rgfilt(:,:,strcmp(dfig1.dname_dmatfilt,dataselected(ii)) & dfig1.rep_dmatfilt == jj);
        tmp = nanmedian(tmp,3);
        medmat4 = cat(3,medmat4,tmp);
    end
    medmats(ii).dmat = medmat;
    medmats(ii).dmatfilt = medmat2(1:end-1,1:end-1,:);
    medmats(ii).Rg = medmat3;
    medmats(ii).Rgfilt = medmat4;
end

%% AVG -- median distancemap 
diagidx = 1:19;
vizidx = [1,3,4,5];%1:length(medmats);%
vizidx = [16,17,18];%1:length(medmats);%
xlabels = {'–75','–30','–15','0','15','30','75'};
xlabelspos = [1,4,7,10,13,16,19];
for ii = 1:length(vizidx)
    fh = SizedFig(20,17);
    medmat = nanmean(medmats(vizidx(ii)).dmatfilt,3);
    imagetriu(rot90(medmat(diagidx,diagidx),2),230,330,flipud(jet),0);
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
figDIR = '/Users/tfuji/Documents/tempanalysis/_figures_v3/fig2/';
diagidx = 1:19;
% vizidx = 1:length(medmats);%[1,3,4,5];
vizidx = [2,3,4,5,7,8,10,12];
vizidx = [17,18];
xlabels = {'–75','–30','–15','0','15','30','75'};
xlabelspos = [1,4,7,10,13,16,19];
% xlabels = {'–30','–15','0','15','30'};
% xlabelspos = [1,4,7,10,13];
for ii = 1:length(vizidx)
    fh = SizedFig(20,17);
%     medmat = nanmean(medmats(vizidx(ii)).dmat - medmats(ctrlidx(vizidx(ii))).dmat,3);
    medmat = nanmean(medmats(vizidx(ii)).dmatfilt - medmats(ctrlidx(vizidx(ii))).dmatfilt,3);
    imagetriu(rot90(medmat(diagidx,diagidx),2),-100,100,flipud(bluewhiteredw0),0);
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
% dataselected = {'KRAB-nodox','KRAB-1day','KRAB-5days','irreversible','reactivated',...
%     'HDAC4-nodox','HDAC4-1day','HDAC4-5days','VP64-nodox','VP64-1day',...
%     'DNMT3B-nodox','DNMT3B-irreversible'};
Nnodes = 5;
ii = 2;
diagidx = 1:13;
medmatd = medmats(ii).Rg(diagidx,diagidx,:);
medmatc = medmats(ctrlidx(ii)).Rg(diagidx,diagidx,:);
localRgD = [];
localRgC = [];
for jj = 1:(size(medmatd,1)- (Nnodes-1))
    localRgD(jj,:) = medmatd(jj,jj+ (Nnodes-1),:);
    localRgC(jj,:) = medmatc(jj,jj+ (Nnodes-1),:);
end
Rgratio = flipud(localRgD./localRgC);
avgRgRatio = nanmean(Rgratio,2);
stdRgRatio = nanstd(Rgratio,[],2) / sqrt(sum(~isnan(Rgratio(1,:))));

fh = SizedFig(30,25);
xax = linspace(-floor(size(localRgD,1)/2)*5,floor(size(localRgD,1)/2)*5,size(localRgD,1));
lh = plot(repmat(xax,[2,1]),[ones(size(localRgD,1),1),avgRgRatio]','b-','linewidth',25);
for ll = 1:length(lh)
    lh(ll).Color=[0,0,1,0.5];
end
hold on;
plot(repmat(xax,[2,1]),[avgRgRatio-stdRgRatio,avgRgRatio+stdRgRatio]','k-','linewidth',1);
% for ii = 1:7
%     plot(xax+(ii-5.5)/5,params(:,ii),'bo','MarkerFaceColor','w','MarkerSize',6);
% end
plot([-30,30],[1.0,1.0],'k-');
plot([0,0],[0.2,2],'k--');
ylim([0.86,1.12]); % for others
% ylim([0.82,1.26]); % for VP64
xlim([-30,30]);
xlabel('position (kb)');
ylabel('Rg ratio (treat/nodox)');
set(gca,'FontSize',15,'XColor','k','YColor','k');
box off;
% saveas(fh,[figDIR,'RgRatio_',dataselected{ii},'.pdf']);


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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lamina distance analysis
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for ii = 1:length(Objs)
    ME = [];
    try
        Objs(ii).df.MaskData;
    catch ME
    end
    if isempty(ME)
        disp(Objs(ii).fname);
    end
end
%%
meddistthresh = 600;
mindist_full = cell(length(Objs),1);
indices = [1,2,4,6,7,8,9];
for kk = indices
    mindist = cell(length(Objs(kk).df),1);

    for j = 1:length(Objs(kk).df)
        % --- need to find indices with repoter DNAFISH signal positive
        tmpcoord = Objs(kk).df(j).coord;
        tmpdmat = Objs(kk).df(j).dmat;
        tmpfilter = repmat(nanmedian(tmpdmat,2) > meddistthresh,[1 3 1]);
        tmpcoord(tmpfilter) = NaN;
        reppos = reshape(~isnan(tmpcoord(10,1,:)),[],1);

        for idx = 1:length( Objs(kk).df(j).MaskData.fiducialx)
            fidx = Objs(kk).df(j).MaskData.fiducialx(idx);
            fidy = Objs(kk).df(j).MaskData.fiducialy(idx);
            edgexy = Objs(kk).df(j).MaskData.edgexy{idx};
            dist = [];
            nan_hybround = isnan(tmpcoord(1:end-1,1,idx));
                if ~isempty(edgexy) && reppos(idx) % && sum(nan_hybround)/19 < 0.5
                    for i = 1:size(edgexy, 1)
                        dist = cat(1, dist, sqrt(sum(([fidx,fidy] - edgexy(i, :)).^2)));
                    end
                end
            mindist{j} = cat(1, mindist{j}, min(dist));
        end
        mindist{j} = mindist{j} * 103;
    end 
    mindist_full{kk} = mindist;
end
%%
conditions = {'KRAB-nodox','KRAB-1day','KRAB-5days','irreversible','reactivated'};
indices = [1,2,4,6,7,8,9];
meddistToEdge = nan(length(conditions),length(Objs));
for mm = 1:length(conditions)
    for jj = 1:length(indices)
        kk = indices(jj);
        for ii = 1:length(Objs(kk).df)
            if strcmp(Objs(kk).df(ii).dname,conditions{mm})
                meddistToEdge(mm,kk) = median(mindist_full{kk}{ii});
            end
        end
    end
end
%% bar graphs for individual condition wiht control
for kk = 2:4
    SizedFig(15,30);
%     subplot(1,3,kk-1);
    hold on;
    tmpmeddist = meddistToEdge([1,kk+1],:)/1000;
    tmpmeddist(1,isnan(tmpmeddist(2,:))) = NaN;
    bar(nanmean(tmpmeddist,2),'FaceColor',[0.8,0.8,0.8]);
    plot(tmpmeddist,'k-');
    for ii = 1:size(tmpmeddist,1)
        tmp = tmpmeddist(ii,:);
        plot(ones(length(tmp),1)*ii,tmp,'bo','MarkerSize',10,'MarkerFaceColor','w',...
            'LineWidth',2);
    end
    title(conditions{kk+1});
    set(gca,'XTick',1:2,'XTickLabel',conditions([1,kk+1]));
    xtickangle(-45);
    ylabel('median distance to the edge of nucleus (um)');
    set(gca,'XColor','k','YColor','k','box','on');
    tmpmeddist = tmpmeddist(:,~isnan(tmpmeddist(1,:)));
    ylim([2,4]);
    [~,pt] = ttest(tmpmeddist(1,:),tmpmeddist(2,:))
    xlim([0.2,2.8]);
end
%%
%%
%%

%%
figure;
indices = [1,2,4,6,7,8,9];
for jj = 1:length(indices)
    hold on;
    kk = indices(jj);
    subplot(3,3,kk);
    for ii = 1:length(Objs(kk).df)
%         if strcmp(Objs(kk).df(ii).dname,'KRAB-nodox') || strcmp(Objs(kk).df(ii).dname,'KRAB-1day')
        if strcmp(Objs(kk).df(ii).dname,'KRAB-nodox') || strcmp(Objs(kk).df(ii).dname,'KRAB-5days')
%         if strcmp(Objs(kk).df(ii).dname,'KRAB-nodox') || strcmp(Objs(kk).df(ii).dname,'irreversible')
%         if strcmp(Objs(kk).df(ii).dname,'KRAB-nodox') || strcmp(Objs(kk).df(ii).dname,'reactivated')
            ecdf(mindist_full{kk}{ii});
        end
    end
    legend(Objs(kk).df.dname);
end

%%
conditions = {'KRAB-nodox','KRAB-1day','KRAB-5days','irreversible','reactivated'};
indices = [1,2,4,6,7,8,9];
meddistToEdge = nan(length(conditions),length(Objs));
for mm = 1:length(conditions)
    for jj = 1:length(indices)
        kk = indices(jj);
        for ii = 1:length(Objs(kk).df)
            if strcmp(Objs(kk).df(ii).dname,conditions{mm})
                meddistToEdge(mm,kk) = median(mindist_full{kk}{ii});
            end
        end
    end
end
%%
figure;
hold on;
bar(nanmean(meddistToEdge,2));
for ii = 1:size(meddistToEdge,1)
    tmp = meddistToEdge(ii,:);
    plot(ones(length(tmp),1)*ii,tmp,'bo','MarkerSize',8,'MarkerFaceColor','w')
end

%%
%%
%%
conditions = {'KRAB-5days','irreversible','reactivated'};
indices = [1,2,4,6,7,8,9];
dist2lamina = [];
dist2lamina.name = [];
dist2lamina.data = [];
for mm = 1:length(conditions)
    cntr = 1;
    for jj = 1:length(indices)
        kk = indices(jj);
        for ii = 1:length(Objs(kk).df)
            if strcmp(Objs(kk).df(ii).dname,conditions{mm})
                idxtmp = Objs(kk).df(ii).ctrlidx;
                dist2lamina(1).data = mindist_full{kk}{idxtmp};
                dist2lamina(1).name = Objs(kk).df(idxtmp).dname;
                dist2lamina(2).data = mindist_full{kk}{ii};
                dist2lamina(2).name = Objs(kk).df(ii).dname;            
                save(['./dist2lamina_',dist2lamina(2).name,'_rep',num2str(cntr),'.mat'],'dist2lamina');
                cntr = cntr + 1;
            end
        end
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% statistical significance test of Rg
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% statistical significance of Rg
close all
SizedFig(15,30);
idx1 = 7-6;
idx2 = 7+6;
hold on;
visidx = [12];%[2,3,4,5,13];%[4];%
disp('----');
for ii = 1:length(visidx)
    tmpidx = reshape(medmats(visidx(ii)).Rg(1,13,:),[],1);
%     tmpidx(3) = NaN;
    tmpidx = ~isnan(tmpidx);
    tmp1 = reshape(medmats(visidx(ii)).Rg(idx1,idx2,:),[],1);
    tmp1 = tmp1(tmpidx);
    tmp2 = reshape(medmats(ctrlidx(visidx(ii))).Rg(idx1,idx2,:),[],1);
    tmp2 = tmp2(tmpidx);
%     xax1 = linspace(ii*2-0.2,ii*2+0.2,length(tmp1));
%     xax2 = linspace(ii*2-1-0.2,ii*2-1+0.2,length(tmp2));
    xax1 = ones(1,length(tmp1))*ii*2;
    xax2 = ones(1,length(tmp1))*(ii*2-1);
    [~,pval] = ttest(tmp1,tmp2);
    disp([num2str(visidx(ii)),':',num2str(pval)]);
    bar(ii*2,nanmean(tmp1),'FaceColor',[0.8,0.8,0.8]);
    bar(ii*2-1,nanmean(tmp2),'FaceColor',[0.8,0.8,0.8]);
    plot([xax1',xax2']',[tmp1,tmp2]','b-','LineWidth',0.5);
    plot(xax1,tmp1,'bo',...
        'MarkerFaceColor','w','MarkerSize',10,'LineWidth',2);
    plot(xax2,tmp2,'bo',...
        'MarkerFaceColor','w','MarkerSize',10,'LineWidth',2);
end
ylabel('med Rg (-30 30) [nm]');
set(gca,'XColor','k','YColor','k','box','on');
ylim([150,220]);
xlim([0.2,2.8]);
%% HDAC4
SizedFig(25,30);
idx1 = 7-6;
idx2 = 7+6;
hold on;
visidx = [1,4,5,12,13];
tmpidx = reshape(medmats(13).Rg(1,13,:),[],1);
tmpidx = ~isnan(tmpidx);
for ii = 1:length(visidx)
    tmp1 = reshape(medmats(visidx(ii)).Rg(1,13,:),[],1);
    tmp1 = tmp1(tmpidx);
    bar(ii,nanmean(tmp1),'FaceColor',[0.8,0.8,0.8]);
%     xax1 = linspace(ii-0.2,ii+0.2,length(tmp1));
    xax1 = ones(1,length(tmp1))*ii;
    plot(xax1,tmp1,'bo',...
        'MarkerFaceColor','w','MarkerSize',10,'LineWidth',2);
    tmp2 = reshape(medmats(ctrlidx(visidx(ii))).Rg(1,13,:),[],1);
    tmp2 = tmp2(tmpidx);
    [~,pval] = ttest(tmp1,tmp2);
    pval
    
    if ii == 4
        tmp2 = reshape(medmats(1).Rg(1,13,:),[],1);
        tmp2 = tmp2(tmpidx);
        [~,pval] = ttest(tmp1,tmp2);
        pval
    end
end
xlim([0.2,length(visidx)+0.8]);
ylabel('med Rg (-30 30) [nm]');
set(gca,'XColor','k','YColor','k','box','on');
ylim([170,210]);ylim([150,220]);
%% HDAC4, with lines
SizedFig(25,30);
idx1 = 7-6;
idx2 = 7+6;
hold on;
visidx = [1,5,4,12,13];
tmpidx = reshape(medmats(13).Rg(1,13,:),[],1);
tmpidx = ~isnan(tmpidx);
for ii = 1:length(visidx)
    tmp1 = reshape(medmats(visidx(ii)).Rg(1,13,:),[],1);
    tmp1 = tmp1(tmpidx);
    bar(ii,nanmean(tmp1),'FaceColor',[0.8,0.8,0.8]);
%     xax1 = linspace(ii-0.2,ii+0.2,length(tmp1));
    xax1 = ones(1,length(tmp1))*ii;
    plot(xax1,tmp1,'bo',...
        'MarkerFaceColor','w','MarkerSize',10,'LineWidth',2);
    tmp2 = reshape(medmats(ctrlidx(visidx(ii))).Rg(1,13,:),[],1);
    tmp2 = tmp2(tmpidx);
    [~,pval] = ttest(tmp1,tmp2);
    pval
    
    if ii == 4
        tmp2 = reshape(medmats(1).Rg(1,13,:),[],1);
        tmp2 = tmp2(tmpidx);
        [~,pval] = ttest(tmp1,tmp2);
        pval
    end
    
    if ii == 2 || ii == 5
        xax2 = ones(1,length(tmp2))*(ii-1);
        plot([xax1',xax2']',[tmp1,tmp2]','b-','LineWidth',0.5);
    end
end
xlim([0.2,length(visidx)+0.8]);
ylabel('med Rg (-30 30) [nm]');
set(gca,'XColor','k','YColor','k','box','on');
ylim([170,210]);ylim([150,220]);
