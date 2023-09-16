%% save data
DIR = '/Users/tfuji/Documents/tempanalysis/_figures_v2/fig2/';
{'KRAB-nodox','KRAB-1day','KRAB-5days','reactivated','irreversible',...
    'KRAB147-nodox','KRAB147-1day','KRAB147-5days','KRAB150-nodox','KRAB150-1day','KRAB150-5days',...
    'HDAC4-nodox','HDAC4-1day','VP64-nodox','VP64-1day'};
ii = 5;
diagidx = 4:16;
medmat = nanmean(medmats(ii).dmatfilt,3);
writematrix(medmat(diagidx,diagidx),[DIR,'dmatfilt_',dataselected{ii},'.csv']);

%% import
DIR = '/Users/tfuji/Documents/tempanalysis/_figures_v2/fig2/supp/randomwalkfit/avg_multiple_study_KRAB-nodox/fittedparams.csv';
paramsN = readmatrix(DIR);
DIR = '/Users/tfuji/Documents/tempanalysis/_figures_v2/fig2/supp/randomwalkfit/avg_multiple_study_KRAB-5days/fittedparams.csv';
paramsK = readmatrix(DIR);
DIR = '/Users/tfuji/Documents/tempanalysis/_figures_v2/fig2/supp/randomwalkfit/avg_multiple_study_irreversible/fittedparams.csv';
paramsI = readmatrix(DIR);
DIR = '/Users/tfuji/Documents/tempanalysis/_figures_v2/fig2/supp/randomwalkfit/avg_multiple_study_reactivated/fittedparams.csv';
paramsR = readmatrix(DIR);

%% bar median, dots individual trials, xreversed
SizedFig(30,30);
hold on;
xax = linspace(-27.5,27.5,12);
params = flipud(paramsN);
plot(xax,median(params,2)*140,'o-','linewidth',1);
params = flipud(paramsK);
plot(xax,median(params,2)*140,'o-','linewidth',1);
ylim([100,160]);
xlim([-30,30]);
xlabel('position (kb)');
ylabel('step size after optimization (nm)');
set(gca,'FontSize',15,'XColor','k','YColor','k');
box on;

%% ------------------------------------------------------------
%% -----------------  random polymer dim. reduction -----------
%% ------------------------------------------------------------
%%
dfall2 = [];
nseg = 13;
npol = 50000;
stepsize = 140;
%
for hh = 1:4
    switch hh 
        case 1
            params = median(paramsN,2);
        case 2
            params = median(paramsK,2);
        case 3
            params = median(paramsR,2);
        case 4
            params = median(paramsI,2);
    end
    rs = stepsize*ones(nseg,1,npol);%normrnd(1,variability,nseg,1,npol);%
    
    for ii = 1:length(params)
        rs(ii+1,:,:) = rs(ii+1,:,:)*params(ii);
    end

    pols = nan(nseg,3,npol);
    thetas = rand(nseg,1,npol)*pi;
    phis   = rand(nseg,1,npol)*2*pi;    
    xs = rs.*sin(thetas).*cos(phis);
    ys = rs.*sin(thetas).*sin(phis);
    zs = rs.*cos(thetas);

    pols(:,1,:) = xs;
    pols(:,2,:) = ys;
    pols(:,3,:) = zs;

    pols = cumsum(pols,1);
    pols = pols - mean(pols,1);
    dfall2(hh).coordinterp = pols;

    % --- visualize dmat
    dmat = zeros(size(pols,1),size(pols,2),size(pols,3));
    for dotpos = 1:size(pols,3)
        tmp = pols(:,:,dotpos);
        for kk = 1:size(tmp,1)
            dmat(:,kk,dotpos) = sqrt(sum((tmp - tmp(kk,:)).^2,2));
        end
    end
    dfall2(hh).dmatinterp = dmat;
end
%% subtracted distance map
% SizedFig(100,20);
SizedFig(100,20);
idxrange = 1:13;
visidx = [1,2,3,4];
xlabels = 1:13;
for ii = 1:length(visidx)
    subplot(1,length(visidx),ii);
    medmat = nanmedian(dfall2(visidx(ii)).dmatinterp,3) - nanmedian(dfall2(1).dmatinterp,3);
    imagetriu(rot90(medmat(idxrange,idxrange),2),-100,100,flipud(bluewhiteredw0),0);
    caxis([-100,100]);
    cbar = colorbar;
    cbar.Label.String = 'nm';
    cbar.Label.FontSize = 15;
    set(gca,'XTick',xlabels,'XTickLabels',xlabels,'FontSize',10);
    title(num2str(visidx(ii)));
end
%% radius of gyration and its fold-change
for ii = 1:length(dfall2)
    tmp = dfall2(ii).coordinterp;
    tmp = tmp - repmat(mean(tmp,1),[size(tmp,1), 1, 1]);
    tmpRg = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
    dfall2(ii).Rg = tmpRg;
end

%% histogram of radius of gyration
cmap = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980];
binnum = 30;
binr = [50,400];
for ii = 1:4
    SizedFig(15,15);
    hold on;
    tmp1 = reshape(dfall2(1).Rg,[],1);
    h1 = histogram(tmp1,binnum,'BinLimits',binr,'normalization','probability',...
        'EdgeAlpha',0.0,'FaceColor',cmap(1,:));
    stairs([h1.BinEdges,h1.BinEdges(end)],[h1.Values,h1.Values(end),0],'k-');

    tmp2 = reshape(dfall2(ii).Rg,[],1);
    h2 = histogram(tmp2,binnum,'BinLimits',binr,'normalization','probability',...
        'EdgeAlpha',0.0,'FaceColor',cmap(2,:));
    stairs([h2.BinEdges,h2.BinEdges(end)],[h2.Values,h2.Values(end),0],'k-');

    plot([median(tmp1),median(tmp1)],[0,0.15],'--','Color',cmap(1,:));
    plot([median(tmp2),median(tmp2)],[0,0.15],'--','Color',cmap(2,:));

    set(gca,'XColor','k','YColor','k');
    median(tmp1)
    median(tmp2)
%     ranksum(tmp1,tmp2)
    box on;
    xlabel('Radius of gyration');
    ylabel('probability');
    ylim([0,0.15]);
end





%% AVG -- median distancemap 
% SizedFig(60,10);
fh = SizedFig(100,20);

diagidx = 4:16;
xlabels = {'–30','–15','0','15','30'};
xlabelspos = [1,4,7,10,13];

vizidx = [1,3,4,5];%1:length(medmats);%

for ii = 1:length(vizidx)
    subplot(1,length(vizidx),ii);
    medmat = nanmean(medmats(vizidx(ii)).dmatfilt,3);
    imagetriu(rot90(medmat(diagidx,diagidx),2),200,350,flipud(jet),0);

    title(dataselected(vizidx(ii)));
    caxis([200,350]);
    cbar = colorbar;
    cbar.Label.String = 'nm';
    cbar.Label.FontSize = 15;
    set(gca,'XTick',xlabelspos,'XTickLabels',xlabels,'FontSize',10,...
        'XColor','k','YColor','k');
    xlabel('kb');
%     saveas(fh,[figDIR,'avg_medmat_',dataselected{vizidx(ii)},'.pdf']);
end




%% AVG -- median distancemap 
% SizedFig(60,10);
fh = SizedFig(100,20);

diagidx = 1:13;
xlabels = {'–30','–15','0','15','30'};
xlabelspos = [1,4,7,10,13];

vizidx = [1,2,3,4];%1:length(medmats);%

for ii = 1:length(vizidx)
    subplot(1,length(vizidx),ii);
    medmat = nanmedian(dfall2(vizidx(ii)).dmatinterp,3);
    imagetriu(rot90(medmat(diagidx,diagidx),2),200,350,flipud(jet),0);

    title(dataselected(vizidx(ii)));
    caxis([200,350]);
    cbar = colorbar;
    cbar.Label.String = 'nm';
    cbar.Label.FontSize = 15;
    set(gca,'XTick',xlabelspos,'XTickLabels',xlabels,'FontSize',10,...
        'XColor','k','YColor','k');
    xlabel('kb');
%     saveas(fh,[figDIR,'avg_medmat_',dataselected{vizidx(ii)},'.pdf']);
end
