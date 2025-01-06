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
% dataselected = {'KRAB-nodox','KRAB-1day','KRAB-5days','reactivated','irreversible',...
%     'KRAB147-nodox','KRAB147-1day','KRAB147-5days','KRAB150-nodox','KRAB150-1day','KRAB150-5days',...
%     'HDAC4-nodox','HDAC4-1day','VP64-nodox','VP64-1day'}; % all
dataselected = {'KRAB-nodox','KRAB-5days','reactivated','irreversible'};%,'KRAB-1day'
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
% dfall = dfall([1:21,23,24]);
% dfall = dfall([1:21,23,24]);
%%
for ii = 1:length(dfall)
    disp([num2str(ii),':',dfall(ii).dname,' -- ',num2str(size(dfall(ii).dmatinterp,3))])
end
%%
for ii = 1:length(dfall)
    disp([num2str(ii),':',dfall(ii).dname,' -- ',num2str(dfall(ii).rep)])
end

%%
pairs_allsample = [];
Rg_allsample = [];
dmat_allsample = [];
dmatnointp_allsample = [];
coord_allsample = [];
coordnointp_allsample = [];
datarange = 1:length(dfall);%
idxrange = 4:16;
for ii = datarange%1:length(dfall)%
    tmpidx = ~isnan(reshape(dfall(ii).coordfilt50p(10,1,:),[],1));
    pairs_allsample = cat(1,pairs_allsample,dfall(ii).pairs);
    Rg_allsample = cat(3,Rg_allsample,dfall(ii).Rg);
    dmat_allsample = cat(3,dmat_allsample,dfall(ii).dmatinterp(idxrange,idxrange,tmpidx));
    dmatnointp_allsample = cat(3,dmatnointp_allsample,dfall(ii).dmatfilt50p(idxrange,idxrange,tmpidx));
    coord_allsample = cat(3,coord_allsample,dfall(ii).coordinterp(idxrange,:,tmpidx));
    coordnointp_allsample = cat(3,coordnointp_allsample,dfall(ii).coordfilt50p(idxrange,:,tmpidx));
end
%% display names of data
for ii = 1:length(dfall)
    disp([num2str(ii),':',dfall(ii).dname]);
end
dataselected = {'KRAB-nodox','KRAB-1day','KRAB-5days','reactivated','irreversible',...
    'KRAB147-nodox','KRAB147-1day','KRAB147-5days','KRAB150-nodox','KRAB150-1day','KRAB150-5days',...
    'HDAC4-nodox','HDAC4-1day','VP64-nodox','VP64-1day'};
ctrlidx = [1,1,1,1,1,6,6,6,9,9,9,12,12,14,14];
IdxByCond = repmat({[]},length(dataselected),1);
for jj = 1:length(dataselected)
    for ii = 1:length(dfall)
        if strcmp(dfall(ii).dname,dataselected{jj})
            IdxByCond{jj} = cat(2,IdxByCond{jj},ii);
        end
    end
end

%% k-means
% rng('default');
kmnum = 4;
% kmnum = 7;
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

%% median of each cluster; x reversed
SizedFig(50,20);
for ii = 1:kmnum
    subplot(1,kmnum,ii);
    tmp = dmatnointp_allsample(:,:,kmopt == ii);
    tmp = rot90(nanmedian(tmp,3),2);
%     tmp = nanmedian(tmp,3);
    imagetriu(tmp,150,330,flipud(jet));
end
%% representative structure of each cluster
    cmapust = [255,0,0; 255,0,0;
        255,131,0; 255,131,0;
        248,255,0; 248,255,0;
    117,255,0; 117,255,0;
    73,255,0; 73,255,0;
    29,255,0; 29,255,0;
    0,255,15; 0,255,15;
    0,255,58; 0,255,58;
    0,255,102; 0,255,102]/255;
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
    
az = 23;%-148;%circshift(1:360,-29);
el = 14;
rotview = 0;
boxwidth = 400;
tmpcoord = coordnointp_allsample;
endpos = ~isnan(tmpcoord(1,1,:)) & ~isnan(tmpcoord(end,1,:));
tmpdmat = dmatnointp_allsample(:,:,endpos);
tmpdmatintp = dmat_allsample(:,:,endpos);
kmopttmp = kmopt(endpos);
tmpcoord = tmpcoord(:,:,endpos);

for jj = 1:kmnum
    SizedFig(20,30);
%     subplot(1,kmnum,jj);
    tmp = tmpdmat(:,:,kmopttmp == jj);
    tmp = nanmedian(tmp,3);
    difffrommed = reshape(nansum(nansum((tmpdmatintp - tmp).^2,1),2),[],1);
    indices = find(difffrommed == min(difffrommed)) % 832,231,1490,2380

    for ii = 1:length(indices)
        tmp = tmpcoord(:,:,indices(ii)); % 616 after thresholding6
            % color code
            cntr = 7;
            ustidx = zeros(size(cmapust,1),1);
            ustidx(1:2:(cntr-1)*2) = flipud(~isnan(tmp(1:(cntr-1),1)));
            ustidx(2:2:(cntr-1)*2) = flipud(~isnan(tmp(1:(cntr-1),1)));
            ustidx = flipud(ustidx);
            dstidx = zeros(size(cmapust,1),1);
            dstidx(1:2:(cntr-1)*2) = ~isnan(tmp((cntr+1):end,1));
            dstidx(2:2:(cntr-1)*2) = ~isnan(tmp((cntr+1):end,1));
            cmap = cat(1,cmapust(ustidx == 1,:),cmaprep,cmapdst(dstidx == 1,:));

        tmp = tmp(~isnan(tmp(:,1)),:); 
        polycenter = (max(tmp,[],1) + min(tmp,[],1))/2;% mean(tmp,1);
        tmp = tmp - polycenter;
 
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

            xbound = [-1,1]*boxwidth;
            ybound = [-1,1]*boxwidth;
            zbound = [-1,1]*boxwidth;
            radius = 30;
            thickness = 10;
    %         clf;%
            hold on;
            [xsphere,ysphere,zsphere] = sphere; 
            xsphere = xsphere*radius; ysphere = ysphere*radius; zsphere = zsphere*radius;
            for i =1:length(tmp)
                % --- plot spheres --- %
                surf(xsphere+tmp(i,1),ysphere+tmp(i,2),zsphere+tmp(i,3),'EdgeAlpha',0,'FaceColor',cmap(2*i-1,:));
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
            axis off;
            box off;
            set(gca,'Projection','perspective','Box','on','BoxStyle','full',...
                'FontSize',10)
            xlim(xbound); ylim(ybound); zlim(zbound);
    %         xlabel('X [nm]');
    %         ylabel('Y [nm]');
    %         zlabel('Z [nm]');
    %         view(az(21),el);
            view(az,el);
    %         view([7.9927 28.5909]);
            camlight('right');
            cntr = 0;
    %         pause(1/10);
        end
    end
end


%% calculate average k_mean percentage for each condition; pi-chart
SizedFig(100,50);
% visidx = [19,23,26,29,30,34];
% % visidx = 1:9;
% % visidx = 10:18;
% visidx = [1,2,3,33,34];
visidx = IdxByCond([1,3,4,5]); % KRABnodox,irrversible,reactivated%
dnum = length(visidx);

kmean_tblall = [];
for ii = 1:dnum
subplot(1,dnum,ii);
tmpidx = visidx{ii};
tmp = [];
    for jj = 1:length(tmpidx)
%         tmp = cat(3,tmp,dfall(tmpidx(jj)).densemat);
        tbl = tabulate(dfall(tmpidx(jj)).km);
        tmp = cat(3,tmp,tbl);
        disp(dfall(tmpidx(jj)).dname);
    end
    pie(mean(tmp(:,3,:),3));
    kmean_tblall = cat(3,kmean_tblall,mean(tmp,3));
    title(dfall(tmpidx(1)).dname);
    hold on;
end

%% k_mean cluster fraction by exp condition and bio rep
visidx = IdxByCond([1,3,4,5]); % IdxByCond(1:15); % KRABnodox,irrversible,reactivated%
dnum = length(visidx);

samplenum = nan(kmnum,dnum,length(Objs));
kmean_tblall = nan(kmnum,dnum,length(Objs));
Rg_all = nan(1,dnum,length(Objs));
for ii = 1:dnum
tmpidx = visidx{ii};
tmp = [];
    for jj = 1:length(tmpidx)
        tbl = tabulate(dfall(tmpidx(jj)).km);
        tmpRg = reshape(nanmedian(dfall(tmpidx(jj)).Rg(1,13,:)),[],1);
        idx = dfall(tmpidx(jj)).rep;
        kmean_tblall(:,ii,idx) = tbl(:,3);
        samplenum(:,ii,idx) = sum(tbl(:,2));
        Rg_all(:,ii,idx) = tmpRg;
    end
end
%% k_mean percentage bar graph / scatter
visidx = [1,2,4,3]; % 1:15;%
conditionnum = length(visidx);
cmap = parula(6);
cmap = cmap(2:end-1,:);
% idxkm = 2; ylims = [0,20];
% idxkm = 3; ylims = [0,20];
idxkm = 1; ylims = [0,70];
% idxkm = 4; ylims = [0,25];
wd = 0.2; 
SizedFig(18,20); 
for i = 1:conditionnum
    kk = visidx(i);
    hold on;
    tmp = reshape(kmean_tblall(idxkm,kk,:),[],1);
    tmpref = reshape(kmean_tblall(idxkm,1,:),[],1);
    tmpref = tmpref(~isnan(tmp));
    tmp = tmp(~isnan(tmp));
    bar(i,nanmean(tmp),0.65,'FaceColor',cmap(idxkm,:),'LineWidth',1);
    xax = linspace(i-wd,i+wd,length(tmp));
    plot(xax,tmp,'ko','MarkerFaceColor','w','MarkerSize',8,...
        'LineWidth',2);
    [~,pv] = ttest(tmpref,tmp);
    disp(pv);
    box on;
end
ylim(ylims);
xlim([0.5,4.5]);
set(gca,'FontSize',20,'XTick',[],'XColor','k','YColor','k');
% xlabel('after dox removal [days]');
ylabel('cluster (%)');
% legend(string(rec_days),'Location','eastoutside');

%% cluster fraction vs median radius of gyration
SizedFig(50,20);
xmin = 150;
xmax = 220;
for ii = 1:4
subplot(1,4,ii);
tmpRg = reshape(Rg_all,[],1);
tmpRg = tmpRg(~isnan(tmpRg));
tmpkm = reshape(kmean_tblall(ii,:,:),[],1);
tmpkm = tmpkm(~isnan(tmpkm));
scatter(tmpRg,tmpkm,80,'filled');
ymin = min(tmpkm)*0.8;
ymax = max(tmpkm)*1.2;
lm = fitlm(tmpRg,tmpkm);
round(lm.Rsquared.Ordinary,2)
intrcpt = lm.Coefficients.Estimate(1);
slope = lm.Coefficients.Estimate(2);
hold on; plot([xmin,xmax],[xmin,xmax]*slope+intrcpt,'k--','LineWidth',0.2);
box on;
set(gca,'XColor','k','YColor','k');
ylim([ymin,ymax]);
xlim([xmin,xmax]);
% corrcoef(tmpRg,tmpkm).^2
end

%% tSNE
% rng('default');
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
SizedFig(40,20);
visidx = [1,2,18,23];
dnum = length(visidx);
for ii = 1:dnum
% subplot(ceil(sqrt(dnum)),ceil(sqrt(dnum)),ii);
subplot(1,dnum,ii);
dscatter(dfall(visidx(ii)).tSNEout(:,1),dfall(visidx(ii)).tSNEout(:,2));
xlabel('t-SNE 1');
ylabel('t-SNE 2');
title(dfall(visidx(ii)).dname);
axis image;
axis([-50,50,-50,50]);
% axis([-60,60,-60,60]);
caxis([0.2*10^-4,0.7*10^-4]);
caxis([0.2*10^-4,0.6*10^-4]);
colormap(jet);

hold on;
box on;
end
%% plot t-sne for individual condition, aggregated reps
SizedFig(40,20);
dataselected = {'KRAB-nodox','KRAB-5days','reactivated','irreversible'};
for ii = 1:length(dataselected)
    % subplot(ceil(sqrt(dnum)),ceil(sqrt(dnum)),ii);
    subplot(1,dnum,ii);
    tsne1_tmp = [];
    tsne2_tmp = [];
    for jj = 1:length(dfall)
        if strcmp(dfall(jj).dname,dataselected{ii})
            tsne1_tmp = [tsne1_tmp;dfall(jj).tSNEout(:,1)];
            tsne2_tmp = [tsne2_tmp;dfall(jj).tSNEout(:,2)];
        end
    end
    dscatter(tsne1_tmp,tsne2_tmp,'MSize',5);
%     histogram(tsne2_tmp);
length(tsne2_tmp)
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    title(dataselected{ii});
    axis image;
    axis([-50,50,-50,50]);
    % axis([-60,60,-60,60]);
    caxis([0.2*10^-4,0.7*10^-4]);
    caxis([0.2*10^-4,0.5*10^-4]);
    colormap(jet);

    hold on;
    box on;
end
%% Tuckey median from bag plot
tmptsne = [];
tmpkm = [];
for ii = 1:length(dfall)%[1,2,3,4,5]%
    tmptsne = cat(1,tmptsne,dfall(ii).tSNEout);
    tmpkm   = cat(1,  tmpkm,dfall(ii).km);
end
clusterMedConv = {};
for ii = 1:4
    tmp = tmptsne(tmpkm == ii,:);
    result = bagplot(tmp,'databag',0,'datafence',0,'plots',0);
    [convidx,~] = convhull(tmp(result.datatype == 2,:));
    tmp = tmp(result.datatype == 2,:);
    clusterMedConv = [clusterMedConv,tmp(convidx,:)];
end
%% plot t-sne for individual condition with Tuckey median
SizedFig(50,30);
visidx = [1,2,18,19];
dnum = length(visidx);
for ii = 1:dnum
subplot(1,dnum,ii);
hold on;
dscatter(dfall(visidx(ii)).tSNEout(:,1),dfall(visidx(ii)).tSNEout(:,2));
plot(clusterMedConv{1}(:,1),clusterMedConv{1}(:,2),'r-','LineWidth',2,'Color','#ff0000');
plot(clusterMedConv{3}(:,1),clusterMedConv{3}(:,2),'r-','LineWidth',2,'Color','#ff22ff');
xlabel('t-SNE 1');
ylabel('t-SNE 2');
title(dfall(visidx(ii)).dname);
axis image;
axis([-50,50,-50,50]);
caxis([0.3*10^-4,0.7*10^-4]);
colormap(parula);

hold on;
box on;
end
%% plot t-sne for individual condition; overlaid
SizedFig(10,15);
idx1 = 1;
idx2 = 3;
tmp1 = dfall(idx1).tSNEout;
tmp2 = dfall(idx2).tSNEout;
tmp = cat(1,tmp1,tmp2);
col1 = [0, 0.4470, 0.7410];
col2 = [0.8500, 0.3250, 0.0980];
cols = [repmat(col1,[size(tmp1,1) 1]);repmat(col2,[size(tmp2,1) 1])];
randidx = randsample(size(tmp,1),size(tmp,1));
scatter(tmp(randidx,1),tmp(randidx,2),6,cols(randidx,:),'filled');
xlabel('t-SNE 1');
ylabel('t-SNE 2');
title(dfall(idx2).dname);
axis image;
axis([-50,50,-50,50]);
set(gca,'XColor','k','YColor','k');

box on;

%% plot t-sne for all condition; colored by radius of gyration
SizedFig(15,20);
tmptsne = [];
for ii = 1:length(dfall)
    tmptsne = cat(1,tmptsne,dfall(ii).tSNEout);
end
scatter(tmptsne(:,1),tmptsne(:,2),5,Rg_allsample,'filled');
xlabel('t-SNE 1');
ylabel('t-SNE 2');
axis image;
axis([-50,50,-50,50]);
caxis([100,250]);
box on;
%% plot t-sne for all condition; colored by k-mean
SizedFig(10,15);
tmptsne = [];
tmpkm = [];
cmap = parula(6);
cmap = cmap(2:end-1,:);
for ii = 1:length(dfall)%[1,2,3,4,5]%
    tmptsne = cat(1,tmptsne,dfall(ii).tSNEout);
    tmpkm   = cat(1,  tmpkm,dfall(ii).km);
end
scatter(tmptsne(:,1),tmptsne(:,2),8,tmpkm,'filled');
colormap(gca,cmap);
xlabel('t-SNE 1');
ylabel('t-SNE 2');
axis image;
axis([-50,50,-50,50]);
set(gca,'XColor','k','YColor','k');

box on;
%% plot t-sne for all condition; colored by k-mean; each cluster separated
SizedFig(40,15);
tmptsne = [];
tmpkm = [];
cmap = parula(6);
cmap = cmap(2:end-1,:);
for ii = 1:length(dfall)%[1,2,3,4,5]%
    tmptsne = cat(1,tmptsne,dfall(ii).tSNEout);
    tmpkm   = cat(1,  tmpkm,dfall(ii).km);
end
for ii = 1:4
    subplot(1,4,ii);
    hold on;
    scatter(tmptsne(tmpkm ~= ii,1),tmptsne(tmpkm ~= ii,2),8,[0.9,0.9,0.9],'filled');
    scatter(tmptsne(tmpkm == ii,1),tmptsne(tmpkm == ii,2),8,cmap(ii,:),'filled');
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    axis image;
    axis([-50,50,-50,50]);
    set(gca,'XColor','k','YColor','k');

    box on;
end
%% average plot for each bin; x reversed
SizedFig(50,20); 
bins = linspace(-50,50,3);
for ii = 1:(length(bins)-1)
    subplot(2,4,ii);
    tmp = dmatnointp_allsample(:,:,tSNEout(:,1) >= bins(ii) & tSNEout(:,1) < bins(ii+1));
    imagetriu(rot90(nanmedian(tmp,3),2),150,350,flipud(jet));
    subplot(2,4,ii+4);
    tmp = dmatnointp_allsample(:,:,tSNEout(:,2) >= bins(ii) & tSNEout(:,2) < bins(ii+1));
    imagetriu(rot90(nanmedian(tmp,3),2),150,350,flipud(jet));
end

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
%%  density 
axrng = [-50,50];
figure;
for idx = 1:length(dfall)
    hh = histogram2(dfall(idx).tSNEout(:,1),dfall(idx).tSNEout(:,2)...
    ,16,'XBinLimits',axrng,'YBinLimits',axrng,'normalization','probability'); 
    dfall(idx).densemat = hh.Values;
    xlabel('x');
    ylabel('y');
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
%     caxis([1.0*10^-5,3.5*10^-5]);
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    hold on;
end

%% average density tsne; subtraction

SizedFig(30,10);
SizedFig(100,10);
% visidx = ([1,2,3,4,5]); % KRABnodox,irrversible,reactivated % IdxByCond;%
visidx = 1:13; % KRABnodox,irrversible,reactivated % IdxByCond;%
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
    caxis([-1,1]*0.4*10^-2);
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    hold on;
    sum(sum(kdematall(:,:,ii).*log(kdematall(:,:,ii)./kdematall(:,:,ctrlidx(ii)))))
end

%% average density tsne; subtraction

SizedFig(40,50);
% visidx = IdxByCond([1,2,3,4,5]); % KRABnodox,irrversible,reactivated
visidx = [1,2,3,6,7,8,9,10,11]; %
dnum = length(visidx);

for ii = 1:dnum
    subplot(ceil(sqrt(dnum)),ceil(sqrt(dnum)),ii);
    tmp = kdematall(:,:,visidx(ii)) - kdematall(:,:,ctrlidx(visidx(ii)));
    imagesc(tmp);
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    title(dataselected(visidx(ii)));
    axis image; axis xy;
    colormap(bluewhiteredw0);
    caxis([-1,1]*0.2*10^-2);
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    hold on;
end
%% radius of gyration
tmp = coord_allsample;
tmp = tmp - repmat(mean(tmp,1),[size(tmp,1), 1, 1]);
Rg_allsample = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
tmpRg = Rg_allsample;
SizedFig(20,40);
tmp1 = tSNEout(:,1);
tmp2 = tSNEout(:,2);
% idx = tmp1 <= 50 & tmp1 >= -50 & tmp2 <= 50 & tmp2 >= -50;
% tmp1 = tmp1(idx);
% tmp2 = tmp2(idx);
% tmpRg = Rg_allsample(idx);
dscatter(tmp2,tmpRg);
corrcoef(tmp2,tmpRg)
set(gca,'YScale','log')
%% Rg histogram by exp condition and bio rep
visidx = IdxByCond;%([1,2,3,4,5]); % KRABnodox,irrversible,reactivated%
dnum = length(visidx);
binnum = 40;
RgHist_all = nan(binnum,dnum,5);
MedRg_all = nan(1,dnum,5);
for ii = 1:dnum
tmpidx = visidx{ii};
tmp = [];
    for jj = 1:length(tmpidx)
        tmp1 = reshape(dfall(tmpidx(jj)).Rg(1,13,:),[],1);
        h1 = histogram(tmp1,binnum,'BinLimits',[50,400],'normalization','probability');
        kltmp = reshape(h1.Values,[],1);
        idx = dfall(tmpidx(jj)).rep;
        RgHist_all(:,ii,idx) = kltmp;
        MedRg_all(:,ii,idx) = median(tmp1);
    end
end
%% KL-divergence by exp condition and bio rep
RgHist_allC = RgHist_all(:,ctrlidx,:);
tmp = RgHist_all.*log(RgHist_all./RgHist_allC);
tmp(isinf(tmp)) = 0;
tmp(isnan(tmp)) = 0;
kldivRg_all = nansum(tmp,1);
kldivRg_all(kldivRg_all == 0) = NaN;
nanmean(kldivRg_all,3)
%% Med Rg Fold-change by exp condition and bio rep
MedRg_allC = MedRg_all(:,ctrlidx,:);
tmp = MedRg_all./MedRg_allC;
tmp(isinf(tmp)) = 0;
tmp(isnan(tmp)) = 0;
MedRgFC_all = nansum(tmp,1);
MedRgFC_all(MedRgFC_all == 0) = NaN;
nanmean(MedRgFC_all,3)

%% histogram of radius of gyration
SizedFig(20,20);
idx1 = 3;
repidx = 2;
tmp1 = RgHist_all(:,ctrlidx(idx1),repidx);
bar(1:length(tmp1),tmp1); alpha(0.5);
hold on;
tmp2 = RgHist_all(:,idx1,repidx);
bar(1:length(tmp2),tmp2); alpha(0.5);
tmp = tmp2.*log(tmp2./tmp1);
tmp(isinf(tmp)) = 0;
tmp(isnan(tmp)) = 0;
nansum(tmp,1)

%% histogram of radius of gyration from specified biorep
cmap = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980];
dnames = {'KRAB5days','reactivated','irreversible'};

for ii = 1:length(dnames)
    fh = SizedFig(15,16);
    hold on;
    switch ii
        case 1 
            idx1 = 1;%19; % krab 5 days
            idx2 = 2;%23; % krab 5 days
        case 2 
            idx1 = 17;%react
            idx2 = 18;%react
        case 3 
            idx1 = 22;%irr
            idx2 = 23;%irr
    end
    binnum = 30;
    binr = [50,400];
    tmp1 = reshape(dfall(idx1).Rg(1,13,:),[],1);
    h1 = histogram(tmp1,binnum,'BinLimits',binr,'normalization','probability',...
        'EdgeAlpha',0.0,'FaceColor',cmap(1,:));
    stairs([h1.BinEdges,h1.BinEdges(end)],[h1.Values,h1.Values(end),0],'k-','LineWidth',.5);

    tmp2 = reshape(dfall(idx2).Rg(1,13,:),[],1);
    h2 = histogram(tmp2,binnum,'BinLimits',binr,'normalization','probability',...
        'EdgeAlpha',0.0,'FaceColor',cmap(2,:));
    stairs([h2.BinEdges,h2.BinEdges(end)],[h2.Values,h2.Values(end),0],'k-','LineWidth',.5);

    plot([median(tmp1),median(tmp1)],[0,0.11],'--','Color',cmap(1,:));
    plot([median(tmp2),median(tmp2)],[0,0.11],'--','Color',cmap(2,:));

    set(gca,'XColor','k','YColor','k');
    median(tmp1)
    median(tmp2)
    ranksum(tmp1,tmp2)
    box on;
    xlabel('Radius of gyration (nm)');
    ylabel('probability');
%     saveas(fh,['./Rg_red=',dnames{ii},'.pdf']);
end


%%
%% k-mean representative structure
tmpkm = [];
tmpcoord = [];
tmpdmat = [];
fig2idx = [IdxByCond{1},IdxByCond{2},IdxByCond{3},IdxByCond{4},IdxByCond{5}];

for ii = 1:length(dfall)%fig2idx%
    tmpkm   = cat(1,  tmpkm,dfall(ii).km);
    tmpcoord   = cat(3,tmpcoord,dfall(ii).coordinterp);
    tmpdmat    = cat(3,tmpdmat ,dfall(ii).dmatinterp);

end
cmap = [117,255,0; 117,255,0;
73,255,0; 73,255,0;
29,255,0; 29,255,0;
0,255,15; 0,255,15;
0,255,58; 0,255,58;
0,255,102; 0,255,102;
220, 165, 104;
0,255,233; 0,255,233;
0,233,255; 0,233,255;
0,189,255; 0,189,255;
0,146,255; 0,146,255;
0,102,255; 0,102,255;
0,58,255; 0,58,255]/255;
for kmidx = 1:4
    SizedFig(35,40); 
%     subplot(1,4,kmidx);
    diagidx = 4:16;
    dmat1 = tmpdmat(diagidx,diagidx,tmpkm == kmidx);
%     dmat1 = tmpdmat(:,:,tmpkm == kmidx);
    coord1 = tmpcoord(:,:,tmpkm == kmidx);
    refmat = nanmedian(dmat1,3);
    % refmat = nanmean(dmat1(:,:,tmptsne(:,2) < 0),3);
    difffrommed = reshape(sum(sum((dmat1 - refmat).^2,1),2),[],1);
    % difffrommed = reshape(sum(sum(abs(dmat1 - refmat),1),2),[],1);
    idxclosest = find(difffrommed == min(difffrommed))
    %
    % figure; 
    % subplot(1,2,1);
    % imagetriu(refmat,200,300,flipud(jet));
    % subplot(1,2,2);
    % imagetriu(dmat1(:,:,idxclosest),0,300,flipud(jet));

    az = 33;
    el = 6;
    xangle = 0; yangle = 0; zangle = 0;
    tmp = coord1(diagidx,:,idxclosest);
    diagidx = 4:16;
    hold on;
    boxhw = 320;    

        tmp = tmp(~isnan(tmp(:,1)),:); 
        % polycenter = mean(tmp,1);
        polycenter = (max(tmp,[],1)+min(tmp,[],1))/2;
        tmp = tmp - polycenter;
        tmp = (rotz(zangle)*roty(yangle)*rotx(xangle)*tmp')';
        tmp(:,3) = tmp(:,3); % to frame the figure in 400x400x400 box. just for aesthetics
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
                        30*200/1000,cmap(i*2,:));
                    % --- optimize axis --- %
                    set(h, 'FaceLighting','phong','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
                    set(h,'EdgeAlpha',0);
                    material shiny; 
                end
            end

        end
    
    % --- optimize axis --- %
        axis equal;
        grid on
        set(gca,'Projection','perspective','Box','off','BoxStyle','full',...
            'FontSize',20)
        xlim(xbound); ylim(ybound); zlim(zbound);
%         xlabel('X [nm]');
%         ylabel('Y [nm]');
%         zlabel('Z [nm]');
        set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],...
            'XColor','k','YColor','k');
    %         view(az(21),el);
        view(az,el);
    %         view([7.9927 28.5909]);
        camlight('right');
        cntr = 0;
        set(gca,'XColor','k','YColor','k','LineWidth',2);%,...
    %         'YGrid','off','XGrid','off','ZGrid','off');
end


%% -------------------------------------------------------------
%% -----------------   emulate   noise     ---------------------
%% -------------------------------------------------------------
%%
fig2idx = [IdxByCond{1},IdxByCond{2},IdxByCond{3},IdxByCond{4},IdxByCond{5}];
%% re-run tSNE
tmptsne = [];
tmpkm = [];
tmpcoord = [];
tmpdmat = [];
tmpcoord50p = [];
tmpallpairs = [];
for ii = 1:length(dfall)%fig2idx%[1,2,3,4,5]%
    tmptsne = cat(1,tmptsne,dfall(ii).tSNEout);
    tmpkm   = cat(1,  tmpkm,dfall(ii).km);
    tmpcoord   = cat(3,tmpcoord,dfall(ii).coordinterp);
    tmpdmat    = cat(3,tmpdmat ,dfall(ii).dmatinterp);
    tmpcoord50p= cat(3,tmpcoord50p,dfall(ii).coordfilt50p);
    tmpallpairs= cat(1,tmpallpairs,dfall(ii).allpairs);
end
%%
% % rng('default')
% tSNEout = tsne(tmpallpairs,'Perplexity',35,'Algorithm','barneshut','Distance','euclidean');

%% re-run tSNE; plot t-sne for individual condition
SizedFig(50,30);

for ii = 1:2

    subplot(1,2,ii);
    switch ii
        case 1
             idx = tmptsne(:,2) < 0;
             dname = 't-SNE low';
        case 2
             idx = tmptsne(:,2) > 0;
             dname = 't-SNE high';
    end
    
    dscatter(tSNEout(idx,1),tSNEout(idx,2));
    xlabel('t-SNE 1'); ylabel('t-SNE 2');
    title(dname);
    axis image; box on;
    axis([-50,50,-50,50]);
    caxis([0.3*10^-4,0.9*10^-4]); colormap(jet); 
    set(gca,'XColor','k','YColor','k');
end

%% re-run tSNE; histogram
cmap = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980];
SizedFig(15,15);
hold on;
binnum = 30;
binr = [-50,50];

for ii = 1:2
    switch ii
        case 1
             idx = tmptsne(:,2) < 0;
             dname = 't-SNE low';
        case 2
             idx = tmptsne(:,2) > 0;
             dname = 't-SNE high';
    end

    tmp1 = tSNEoutN(idx,2);
    h1 = histogram(tmp1,binnum,'BinLimits',binr,'normalization','probability',...
        'EdgeAlpha',0.0,'FaceColor',cmap(ii,:));
    stairs([h1.BinEdges,h1.BinEdges(end)],[h1.Values,h1.Values(end),0],'k-');

    set(gca,'XColor','k','YColor','k');

    xlabel('t-SNE 1');
    ylabel('probability');
    box on;
end


%% generate probability distribution
dataselected = {'KRAB-nodox','KRAB-1day','KRAB-5days','reactivated','irreversible',...
    'KRAB147-nodox','KRAB147-1day','KRAB147-5days','KRAB150-nodox','KRAB150-1day','KRAB150-5days',...
    'HDAC4-nodox','HDAC4-1day','HDAC4-5days','VP64-nodox','VP64-1day','DNMT3B-nodox','DNMT3B-irreversible'}; % all
dfallcoordfilt = [];
for idxobj = 1:length(Objs) % [1,2,3,4,5,6,7,10] % 
    for jj = 1:length(Objs(idxobj).df)
%          if sum(strcmp(Objs(idxobj).df(jj).dname,dataselected))
            try
                tmp = rmfield(Objs(idxobj).df(jj),'connmap');
            catch ME
                tmp = Objs(idxobj).df(jj);
            end
            dfallcoordfilt = cat(2,dfallcoordfilt,tmp);
%          end
    end
end
rehyb = [];
for idx = 1:length(dfallcoordfilt)%[1,2,3,4,5]%1:9
    tmp = reshape(dfallcoordfilt(idx).coordfilt(7,:,:) - dfallcoordfilt(idx).coordfilt(20,:,:),3,[])';
    rehyb = cat(1,rehyb,tmp);
end
rehyb = rehyb(~isnan(rehyb(:,1)),:);
SizedFig(30,10);
pd = {};
pd2 = {};
xtitles = {'Xerror (nm)','Yerror (nm)','Zerror (nm)'};
for ii = 1:3
    subplot(1,3,ii);
    pd{ii} = fitdist(rehyb(:,ii),'tLocationScale');
%     pd2{ii}= makedist('tLocationScale','mu',0.0,'sigma',pd{ii}.sigma,'nu',pd{ii}.nu);
    hh = histogram(rehyb(:,ii),100,'BinLimits',[-800,800],'normalization','probability','EdgeColor','none');
    hold on;
    x_values = linspace(-800,800,100);
    y = pdf(pd{ii},x_values);
    y = y/sum(y);
    xlabel(xtitles{ii});
    pd{ii}.sigma
    plot(x_values,y,'LineWidth',1);
    set(gca,'XColor','k','YColor','k');
end
%% add noise
% rng('default') % For reproducibility
tmpcoordN = tmpcoord;
tmpcoordN50p = tmpcoord50p;
for ii = 1:3
    tmp = tmpcoord(:,ii,:);
    rn = random(pd{ii},size(tmp));
    tmpcoordN(:,ii,:) = tmp + rn;
    tmp = tmpcoord50p(:,ii,:);
    rn = random(pd{ii},size(tmp));
    tmpcoordN50p(:,ii,:) = tmp + rn;
end
tmpd = tmpcoordN50p(1:end-1,:,:);
tmpnand = nan(31,3,size(tmpd,3));
tmpnand([1,4,7,10:22,25,28,31],:,:) = tmpd;
coordinterp = fillmissing(tmpnand,'linear',1);
coordinterp = coordinterp([1,4,7,10:22,25,28,31],:,:);
dmatinterp = zeros(size(coordinterp,1),size(coordinterp,2),size(coordinterp,3));
for dotpos = 1:size(coordinterp,3)
    tmp = coordinterp(:,:,dotpos);
    for kk = 1:size(tmp,1)
        dmatinterp(:,kk,dotpos) = sqrt(sum((tmp - tmp(kk,:)).^2,2));
    end
end
% calculate dmat
distancemat = zeros(size(tmpcoordN,1),size(tmpcoordN,1),size(tmpcoordN,3));
for dotpos = 1:size(tmpcoordN,3)
    tmp2 = tmpcoordN(:,:,dotpos);
    for kk = 1:size(tmp2,1)
        distancemat(:,kk,dotpos) = sqrt(sum((tmp2 - tmp2(kk,:)).^2,2));
    end
end

%% visualize 3d structure
boxwidth = 800;
% idxseries = [832,231,1490,2380];
% idxseries = [1097,177,1258,1936];
idxseries = [364,1184,226,1475];
diagidx = 4:16;
repcoord = tmpcoord50p(diagidx,:,:);
repcoordN = tmpd(diagidx,:,:);
endpos = ~isnan(repcoord(1,1,:)) & ~isnan(repcoord(end,1,:));
repcoord = repcoord(:,:,endpos);
repcoordN = repcoordN(:,:,endpos);

for kmidx = 1:4
    SizedFig(20,30); 


    az = 23;%-148;%circshift(1:360,-29);
    el = 14;
    hold on;
    for ii = 1:2
        switch ii
            case 1
                tmp = repcoord(:,:,idxseries(kmidx)); 
                cmap = [0, 0.4470, 0.7410];
                tmp = tmp(~isnan(tmp(:,1)),:); 
                polycenter = (max(tmp,[],1)+min(tmp,[],1))/2;
                tmp = tmp - polycenter;
            case 2
                tmp = repcoordN(:,:,idxseries(kmidx));  
                cmap = [0.8500, 0.3250, 0.0980];
                tmp = tmp(~isnan(tmp(:,1)),:); 
                tmp = tmp - polycenter;
        end




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

            xbound = [-1,1]*boxwidth;
            ybound = [-1,1]*boxwidth;
            zbound = [-1,1]*boxwidth;
            radius = 30;
            thickness = 10;
            hold on;
            [xsphere,ysphere,zsphere] = sphere; 
            xsphere = xsphere*radius; ysphere = ysphere*radius; zsphere = zsphere*radius;
            for i =1:length(tmp)
                % --- plot spheres --- %
                surf(xsphere+tmp(i,1),ysphere+tmp(i,2),zsphere+tmp(i,3),'EdgeAlpha',0,'FaceColor',cmap);
                if i ~= length(tmp)
                    % --- plot spline interpolation --- %
                    pltrange = (1+10*(i-1)):(1+10*i);
                    h = plot3t(tmpinterpx(pltrange),tmpinterpy(pltrange),tmpinterpz(pltrange),...
                        thickness,cmap);
                    % --- optimize axis --- %
                    set(h, 'FaceLighting','phong','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
                    set(h,'EdgeAlpha',0);
                    material shiny; 
                end
            end

        end
    end
    % --- optimize axis --- %
            axis equal;
            grid off;
            axis off;
            box off;
            set(gca,'Projection','perspective','Box','on','BoxStyle','full',...
                'FontSize',10)
            xlim(xbound); ylim(ybound); zlim(zbound);
    %         xlabel('X [nm]');
    %         ylabel('Y [nm]');
    %         zlabel('Z [nm]');
    %         view(az(21),el);
            view(az,el);
    %         view([7.9927 28.5909]);
            camlight('right');
            cntr = 0;
    %         pause(1/10);
end
%% tSNE; histogram Rg
diagidx = 4:16;
% diagidx = 1:13;
repcoord = tmpcoord(diagidx,:,:);
Rg1 = reshape(sqrt(mean(sum((repcoord - mean(repcoord,1)).^2,2))),[],1);
repcoordN = coordinterp(diagidx,:,:);
Rg2 = reshape(sqrt(mean(sum((repcoordN - mean(repcoordN,1)).^2,2))),[],1);
% endpos = ~isnan(repcoord(1,1,:)) & ~isnan(repcoord(end,1,:));
% repcoord = repcoord(:,:,endpos);
% repcoordN = repcoordN(:,:,endpos);

cmap = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980];
binnum = 30;
binr = [0,500];
visidx = [1,3];
for ii = 1:2
SizedFig(15,15);
hold on;
    switch ii 
        case 1
            tmp1 = Rg1(tmptsne(:,2) > 0);
            tmp1 = tmp1(randsample(size(tmp1,1),1363));
            tmp2 = Rg1(tmptsne(:,2) < 0);
            tmp2 = tmp2(randsample(size(tmp2,1),2566));
        case 2
            tmp1 = Rg2(tmptsne(:,2) > 0);
            tmp1 = tmp1(randsample(size(tmp1,1),1363));
            tmp2 = Rg2(tmptsne(:,2) < 0);
            tmp2 = tmp2(randsample(size(tmp2,1),2566));
    end
    
    h1 = histogram(tmp1,binnum,'BinLimits',binr,'normalization','probability',...
        'EdgeAlpha',0.0,'FaceColor',cmap(ii,:)*0.5);
    stairs([h1.BinEdges,h1.BinEdges(end)],[h1.Values,h1.Values(end),0],'k-');

    h1 = histogram(tmp2,binnum,'BinLimits',binr,'normalization','probability',...
        'EdgeAlpha',0.0,'FaceColor',cmap(ii,:));
    stairs([h1.BinEdges,h1.BinEdges(end)],[h1.Values,h1.Values(end),0],'k-');

    set(gca,'XColor','k','YColor','k');
    ranksum(tmp1,tmp2)
    xlabel('radius of gyration (nm)');
    ylabel('probability');
    box on;
end
%% calc all pair dist after adding noise
tmpallpairsN = [];
for ii = 1:size(coordinterp,3)
    tmp = coordinterp(4:16,:,ii);
    tmpd = [];
    for jj = 1:(size(tmp,1)-1)
        for kk = (jj+1):size(tmp,1)
            tmpd = [tmpd,sqrt(sum((tmp(jj,:) - tmp(kk,:)).^2,2))];
        end
    end
    tmpallpairsN = [tmpallpairsN;tmpd];
end

%% run t-SNE after adding noise
rng('default')
tSNEoutN = tsne(tmpallpairsN,'Perplexity',35,'Algorithm','barneshut','Distance','euclidean');

%% tSNE; scatter plot
cmap = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980];

for ii = 1:2
SizedFig(15,15);
hold on;
    switch ii 
        case 1
            tmp1 = tmptsne(tmptsne(:,2) > 0,:);
            tmp2 = tmptsne(tmptsne(:,2) < 0,:);
        case 2
            tmp1 = tSNEoutN(tmptsne(:,2) > 0,:);
            tmp2 = tSNEoutN(tmptsne(:,2) < 0,:);
    end
    subplot(1,2,1);
    plot(tmp1(:,1),tmp1(:,2),'.','Color',cmap(ii,:)*0.5);
    subplot(1,2,2);
    plot(tmp2(:,1),tmp2(:,2),'.','Color',cmap(ii,:));

    set(gca,'XColor','k','YColor','k');
%     xlabel('radius of gyration (nm)');
%     ylabel('probability');
    box on;
end

%% tSNE; scatter plot
cmap = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980];

for ii = 1:2
SizedFig(15,15);
hold on;
    switch ii 
        case 1
            tmp1 = tmptsne(tmptsne(:,2) > 0,:);
            tmp2 = tmptsne(tmptsne(:,2) < 0,:);
        case 2
            tmp1 = tSNEoutN(tmptsne(:,2) > 0,:);
            tmp2 = tSNEoutN(tmptsne(:,2) < 0,:);
    end
    subplot(1,2,1);
    dscatter(tmp1(:,1),tmp1(:,2));
    box on; axis image;
    axis([-50,50,-50,50]);

    subplot(1,2,2);
    dscatter(tmp2(:,1),tmp2(:,2));
    set(gca,'XColor','k','YColor','k');
%     xlabel('radius of gyration (nm)');
%     ylabel('probability');
    box on; axis image;
    axis([-50,50,-50,50]);
end
%% AVG -- median distancemap 
% export medmat from fig1.m 
fh = SizedFig(25,20);


diagidx = 1:20;

vizidx = [1];%1:length(medmats);%
xlabels = [1,4,7,10,13,16,19];
xlabelspos = [1,4,7,10,13,16,19];
for ii = 1:length(vizidx)
    subplot(1,length(vizidx),ii);
    medmat = nanmean(medmats(vizidx(ii)).dmatfilt(:,:,1),3);
    imagetriu(medmat(diagidx,diagidx),100,300,flipud(jet),0);
    title(dataselected(vizidx(ii)));
    caxis([100,300]);
    cbar = colorbar;
    cbar.Label.String = 'nm';
    cbar.Label.FontSize = 15;
    set(gca,'XTick',xlabelspos,'XTickLabels',xlabels,'FontSize',10,...
        'XColor','k','YColor','k');
    xlabel('hyb round');
%     saveas(fh,[figDIR,'avg_medmat_',dataselected{vizidx(ii)},'.pdf']);
end

%%
%%
nurange = [0.5,1,1.5,2,4,8,10,12];
sigrange = [10,20,30,40,80,120,160,200];
sigmas = zeros(length(nurange),length(sigrange),2);
nus    = zeros(length(nurange),length(sigrange),2);
for ii = 1:length(nurange)
    for jj = 1:length(sigrange)
        pdtmp = makedist('tLocationScale','mu',0.0,'sigma',sigrange(jj),'nu',nurange(ii));
        rn1 = random(pdtmp,400000,1);
        rn2 = random(pdtmp,400000,1);
        pdtmp1 = fitdist(rn1,'tLocationScale');
        pdtmp2 = fitdist(rn1-rn2,'tLocationScale');
        nus(ii,jj,1) = pdtmp1.nu;
        nus(ii,jj,2) = pdtmp2.nu;
        sigmas(ii,jj,1) = pdtmp1.sigma;
        sigmas(ii,jj,2) = pdtmp2.sigma;
    end
end
%%
figure;
subplot(1,2,1);
plot(nus(:,:,1)',nus(:,:,2)','o');
subplot(1,2,2);
plot(sigmas(:,:,1)',sigmas(:,:,2)','o');
hold on;
slopes =[];
for ii = 1:8
    lm = fitlm(sigmas(ii,:,1)',sigmas(ii,:,2)');
    intcpt = lm.Coefficients.Estimate(1);
    slope = lm.Coefficients.Estimate(2);
    plot([min(sigmas(ii,:,1)'),max(sigmas(ii,:,1)')],...
        [min(sigmas(ii,:,1)'),max(sigmas(ii,:,1)')]*slope+intcpt,'-');
    slopes = [slopes;slope];
end
% nu2 = nu*1.73
% sigma2 = 
%%
fitlm(reshape(nus(:,:,1),[],1),reshape(nus(:,:,2),[],1))

%%
%%
%%
%% rank sum for each dataset
cmap = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980];
dnames = {'KRAB-5days','reactivated','irreversible'};
rksm = [];
for kk = 1
    for ii = 1:length(dfall)
        if strcmp(dfall(ii).dname, dnames{kk})
%             print('1');
            tmp1 = reshape(dfall(ii).Rg(1,13,:),[],1);
            for jj = 1:length(dfall)
                if strcmp(dfall(jj).dname, 'KRAB-nodox') && dfall(jj).rep == dfall(ii).rep
                    tmp2 = reshape(dfall(jj).Rg(1,13,:),[],1);
                end
            end
            disp([length(tmp1),length(tmp2)]);
            rksm = [rksm;ranksum(tmp1,tmp2)];
        end
    end
end 



%% plot Rg for individual condition, aggregated reps
SizedFig(40,20);
dataselected = {'KRAB-nodox','KRAB-5days','reactivated','irreversible'};
for ii = 2%:length(dataselected)
    % subplot(ceil(sqrt(dnum)),ceil(sqrt(dnum)),ii);
    subplot(1,dnum-1,ii-1);
    Rg_tmp = [];
    Rg_tmp_ctrl = [];
    for jj = 1:2% 1:length(dfall)
        if strcmp(dfall(jj).dname,dataselected{ii})
            Rg_tmp = [Rg_tmp;reshape(dfall(jj).Rg(1,13,:),[],1)];
        end
    end
    for jj = 1:2% 1:length(dfall)
        if strcmp(dfall(jj).dname,dataselected{1})
            Rg_tmp_ctrl = [Rg_tmp_ctrl;reshape(dfall(jj).Rg(1,13,:),[],1)];
        end
    end
    
    histogram(Rg_tmp_ctrl,40,'BinLimits',[0,700],'normalization','probability');
    hold on;
    histogram(Rg_tmp,40,'BinLimits',[0,700],'normalization','probability');
%     ecdf(Rg_tmp_ctrl);
%     hold on;
%     ecdf(Rg_tmp);
    disp([size(Rg_tmp,1),size(Rg_tmp_ctrl,1)]);
    ranksum(Rg_tmp_ctrl,Rg_tmp)
%     histogram(tsne2_tmp);
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    title(dataselected{ii});
    xlim([0,500]);
    % axis([-60,60,-60,60]);
    hold on;
    box on;
end
%%
%%
rksm_K5 = log10(1.5728e-12);
ranksum_reps = [];
for ii = 1:100
    tmp1 = Rg2(tmptsne(:,2) > 0);
    tmp1 = tmp1(randsample(size(tmp1,1),1363));
    tmp2 = Rg2(tmptsne(:,2) < 0);
    tmp2 = tmp2(randsample(size(tmp2,1),2566));
    ranksum_reps = [ranksum_reps;ranksum(tmp1,tmp2)];
end

figure;
hold on;
plot(0,rksm_K5,'bo');
plot([-0.1,0.1],[log10(ranksum_reps),log10(ranksum_reps)],'r-');
plot([-0.1,0.1],[1,1],'k--');
















