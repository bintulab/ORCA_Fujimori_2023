
%% directory for figure
figDIR = '/Users/tfuji/Documents/tempanalysis/_figures_v3/fig4_5/';
%% re-aggregate data

dataselected = {'KRAB-nodox','KRAB-1day','KRAB-5days','irreversible','reactivated',...
    'KRAB147-nodox','KRAB147-1day','KRAB147-5days',...
    'HDAC4-nodox','HDAC4-1day','HDAC4-5days','DNMT3B-nodox','DNMT3B-irreversible',...
    'KRAB150-nodox','KRAB150-1day','KRAB150-5days','VP64-nodox','VP64-1day'};
ctrlidx = [1,1,1,1,1,6,6,6,9,9,9,12,12,14,14,14,17,17];
%%
dfig2 = [];
dfig2.dmat = [];
dfig2.dmatNoIntp = [];
dfig2.dmatfilt = [];
dfig2.coordfilt = [];
dfig2.coordNoIntp = [];
dfig2.Rg = [];
dfig2.Rgfilt = [];
dfig2.rep = [];
dfig2.rep_dmatfilt = [];
dfig2.ctrlidx = [];
dfig2.ctrlidx_dmatfilt = [];
dfig2.dname = {};
dfig2.dname_dmatfilt = {};
for kk = 1:length(Objs)
   df1 = Objs(kk).df;
    drange = 1:length(df1);
    for ii = 1:length(drange)
        if sum(strcmp(df1(ii).dname,dataselected))
            dnum = size(df1(ii).dmatinterp,3);
            dfig2.dmat = cat(3,dfig2.dmat,df1(ii).dmatinterp);
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
            dfig2.Rg = cat(3,dfig2.Rg,tmpRg);
            dfig2.dmatNoIntp = cat(3,dfig2.dmatNoIntp,df1(ii).dmatfilt50p);
            dfig2.coordNoIntp = cat(3,dfig2.coordNoIntp,df1(ii).coordfilt50p);
            dfig2.dmatfilt = cat(3,dfig2.dmatfilt,df1(ii).dmatfilt);
            dfig2.coordfilt = cat(3,dfig2.coordfilt,df1(ii).coordfilt);
            for jj = 1:dnum
                dfig2.rep = cat(1,dfig2.rep,kk);            
                dfig2.ctrlidx = cat(1,dfig2.ctrlidx,ctrlidx(strcmp(df1(ii).dname,dataselected)));            
                dfig2.dname = [dfig2.dname;df1(ii).dname];            
            end
            Objs(kk).df(ii).Rg = tmpRg;
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
            dfig2.Rgfilt = cat(3,dfig2.Rgfilt,tmpRg);
            for jj = 1:size(df1(ii).dmatfilt,3)
                dfig2.dname_dmatfilt = [dfig2.dname_dmatfilt;df1(ii).dname];            
                dfig2.rep_dmatfilt = cat(1,dfig2.rep_dmatfilt,kk);            
            end
        end
    end
end


%% -----------------------------------------------------------------
%% --------- calculate median distance map for each bioreps --------
%% -----------------------------------------------------------------
%% re-group data by experimental condition
% SizedFig(40,10);
c_thresh = 200; 
diagidx = 4:16;
% diagidx = 1:13;
% diagidx = 1:19;
xlabels = {'–30','–20','–10','0','10','20','30'};
medmats = [];
medmats.dmat = [];
medmats.dmatNoIntp = [];
for ii = 1:length(dataselected)
    medmat = [];
    medmat2 = [];
    medmat3 = [];
    medmat4 = [];
    medmat6 = [];
    for jj = 1:length(Objs)
        tmp = dfig2.dmat(:,:,strcmp(dfig2.dname,dataselected(ii)) & dfig2.rep == jj);
        tmp = nanmedian(tmp,3);
        medmat = cat(3,medmat,tmp);
        tmp = dfig2.dmatfilt(:,:,strcmp(dfig2.dname_dmatfilt,dataselected(ii)) & dfig2.rep_dmatfilt == jj);
        tmp = nanmedian(tmp,3);
        medmat2 = cat(3,medmat2,tmp);
        tmp = dfig2.Rg(:,:,strcmp(dfig2.dname,dataselected(ii)) & dfig2.rep == jj);
        tmp = nanmedian(tmp,3);
        medmat3 = cat(3,medmat3,tmp);
        tmp = dfig2.Rgfilt(:,:,strcmp(dfig2.dname_dmatfilt,dataselected(ii)) & dfig2.rep_dmatfilt == jj);
        tmp = nanmedian(tmp,3);
        medmat4 = cat(3,medmat4,tmp);
        tmp = dfig2.dmatfilt(:,:,strcmp(dfig2.dname_dmatfilt,dataselected(ii)) & dfig2.rep_dmatfilt == jj);
        tmpcfreq = nansum(tmp < c_thresh,3)./sum(~isnan(tmp),3);
        medmat6 = cat(3,medmat6,tmpcfreq);
    end
    medmats(ii).dmat = medmat;
    medmats(ii).dmatfilt = medmat2(1:end-1,1:end-1,:);
    medmats(ii).Rg = medmat3;
    medmats(ii).Rgfilt = medmat4;
    medmats(ii).cfreq = medmat6;
end
%
Rg_series = [];
idx_st = 1; idx_ed = 13;
for ii = 1:length(medmats)
    Rg_series = [Rg_series,reshape(medmats(ii).Rg(idx_st,idx_ed,:),[],1)];
end

%% AVG -- median distancemap 
% SizedFig(60,10);
fh = SizedFig(100,20);

diagidx = 1:19;
diagidx = 4:16;
% diagidx = 1:13;
% diagidx = [1:3,4,7,10,13,16,17:19];

vizidx = [2,3,7,8];%1:length(medmats);%
% vizidx = [1,2,3,7,8];%1:length(medmats);%
% vizidx = [9,10,11];%1:length(medmats);%
% xlabels = {'–30','–20','–10','0','10','20','30'};
xlabels = {'–75','–30','–15','0','15','30','75'};
xlabelspos = [1,4,7,10,13,16,19];
for ii = 1:length(vizidx)
    subplot(1,length(vizidx),ii);
    medmat = nanmean(medmats(vizidx(ii)).cfreq,3);
%     imagetriu(rot90(medmat(diagidx,diagidx),2),230,330,flipud(jet),0);
    imagetriu(rot90(medmat(diagidx,diagidx),2),0,0.5,jet,0);
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


% diagidx = [1,2,3,4,7,10,13,16,17,18,19];%1:19;
diagidx = 1:13;
% diagidx = 4:16;
% diagidx = 1:19;

% vizidx = [2,4,6];
% vizidx = [7,8,10,11,13];
% vizidx = 1:length(medmats);%[1,3,4,5];
vizidx = [2,3,7,8,15,16];
% vizidx = 18;
xlabels = {'–75','–30','–15','0','15','30','75'};
xlabelspos = [1,4,7,10,13,16,19];
% xlabels = {'–30','–15','0','15','30'};
% xlabelspos = [1,4,7,10,13];
for ii = 1:length(vizidx)
    fh = SizedFig(20,15);
%     subplot(1,length(vizidx),ii);
    medmat = nanmean(medmats(vizidx(ii)).Rg - medmats(ctrlidx(vizidx(ii))).Rg,3);
%     medmat = nanmean(medmats(vizidx(ii)).dmat - medmats(1).dmat,3);
%     medmat = nanmean(medmats(vizidx(ii)).dmatfilt - medmats(ctrlidx(vizidx(ii))).dmatfilt,3);
%     medmat = nanmean(log(medmats(vizidx(ii)).dmatfilt./medmats(ctrlidx(vizidx(ii))).dmatfilt),3);
%     for kk = 1:size(medmat,1)
%         medmat(kk,kk) = 0;
%     end
%     medmat = nanmean(medmats(vizidx(ii)).dmatfilt - medmats(1).dmatfilt,3);
%     medmat(medmat > -8 & medmat < 5) = 0;
    imagetriu(rot90(medmat(diagidx,diagidx),2),-30,30,flipud(parula),0);
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

%% contact frequency
fh = SizedFig(60,13); 
hold on;
diagidx = 4:16;
% diagidx = 1:19;
ctrlcol = [0.6,0.6,0.6];
vizidx = [7,8,2,3];
% vizidx = [1,6,15,16]; % control
for ii = 1:length(vizidx)
    subplot(1,4,ii);
    hold on;
    tmp1 = flipud(medmats(vizidx(ii)).cfreq(diagidx,10,:));
    tmp2 = flipud(medmats(ctrlidx(vizidx(ii))).cfreq(diagidx,10,:));
    avgtmp1 = nanmean(tmp1,3);
    avgtmp1 = mean(avgtmp1([1:6,8:13]));
    avgtmp2 = nanmean(tmp2,3);
    avgtmp2 = mean(avgtmp2([1:6,8:13]));
    
    plot(nanmean(tmp2,3),'ro-','Color',ctrlcol)
    plot([1:length(diagidx);1:length(diagidx)],squeeze(tmp2(:,:,1:2))','-','Color',ctrlcol)
    plot([1,length(diagidx)],[avgtmp2,avgtmp2],'r--','Color',ctrlcol)

    plot(nanmean(tmp1,3),'bo-')
    plot([1:length(diagidx);1:length(diagidx)],squeeze(tmp1(:,:,1:2))','b-')
    plot([1,length(diagidx)],[avgtmp1,avgtmp1],'b--')
    ylim([0.2,0.59]);
    xlim([0.5,13.5]);
    set(gca,'XTick',[1,4,7,10,13],'XTickLabel',[-30,-15,0,15,30],'XColor','k','YColor','k');
%     set(gca,'XTick',[1,4,7,10,13],'XColor','k','YColor','k');
    xlabel('kb');
    ylabel('contact frequency');
    title(dataselected{vizidx(ii)});
end
% saveas(fh,[figDIR,'cfreq_memoryseries_.pdf']);
% saveas(fh,[figDIR,'cfreq_memoryseries_ctrl.pdf']);
%% contact frequency with power law fit

xax = (0:5:30)';
xax_o = (0:0.1:30)';
% Assuming your data is in arrays: distance and contact_freq
ft = fittype('a*abs(x)^(-b)');
% Create fit options
opts = fitoptions(ft);
opts.StartPoint = [1, 1];  % Initial guess for parameters
opts.Lower = [1, 0];         % Lower bounds
opts.Upper = [Inf, Inf];     % Upper bounds
plawfit_b = zeros(8,2);
xaxp = 1:13;
fh = SizedFig(60,13); 
hold on;
diagidx = 4:16;
% diagidx = 1:19;
ctrlcol = [0.6,0.6,0.6];
% vizidx = [7,8,2,3];
vizidx = [4,5,10,11,13,18];
% vizidx = [1,6,15,16]; % control
for ii = 1:length(vizidx)
    subplot(1,length(vizidx),ii);
    hold on;
    tmp1 = flipud(medmats(vizidx(ii)).cfreq(diagidx,10,:));
    tmp2 = flipud(medmats(ctrlidx(vizidx(ii))).cfreq(diagidx,10,:));
    avgtmp1 = nanmean(tmp1,3);
    avgtmp1 = mean(avgtmp1([1:6,8:13]));
    avgtmp2 = nanmean(tmp2,3);
    avgtmp2 = mean(avgtmp2([1:6,8:13]));
    
    plot(nanmean(tmp2,3),'ro','Color',ctrlcol)
%     plot([1:length(diagidx);1:length(diagidx)],squeeze(tmp2(:,:,1:2))','-','Color',ctrlcol)
    plot([1:length(diagidx);1:length(diagidx)],[nanmean(tmp2,3)'-nanstd(tmp2,[],3)'; nanmean(tmp2,3)'+nanstd(tmp2,[],3)'],'-','Color',ctrlcol);
    
%     plot([1,length(diagidx)],[avgtmp2,avgtmp2],'r--','Color',ctrlcol)
        % cfreq fit
        cfreq_prof = nanmean(tmp2,3);
        contact_freq = cfreq_prof(7:end);
        [fitted_curve_r, goodness] = fit(xax(2:end), contact_freq(2:end), ft, opts);
        contact_freq = flipud(cfreq_prof(1:7));
        [fitted_curve_l, goodness] = fit(-xax(2:end), contact_freq(2:end), ft, opts);
        plot(xax_o/30*(13-7)+7,fitted_curve_r(xax_o),'-','Color',ctrlcol)
        plot(-xax_o/30*(13-7)+7,fitted_curve_l(xax_o),'-','Color',ctrlcol)
        disp('ctrl');
        fprintf('Decay constant: %.3f\n', fitted_curve_r.b)
        fprintf('Decay constant: %.3f\n', fitted_curve_l.b)
        plawfit_b(ii*2-1,1) = fitted_curve_l.b;
        plawfit_b(ii*2-1,2) = fitted_curve_r.b;
    
    plot(nanmean(tmp1,3),'bo')
%     plot([1:length(diagidx);1:length(diagidx)],squeeze(tmp1(:,:,1:2))','b-')
    plot([1:length(diagidx);1:length(diagidx)],[nanmean(tmp1,3)'-nanstd(tmp1,[],3)'; nanmean(tmp1,3)'+nanstd(tmp1,[],3)'],'b-');
%     plot([1,length(diagidx)],[avgtmp1,avgtmp1],'b--')
        % cfreq fit
        cfreq_prof = nanmean(tmp1,3);
        contact_freq = cfreq_prof(7:end);
        [fitted_curve_r, goodness] = fit(xax(2:end), contact_freq(2:end), ft, opts);
        contact_freq = flipud(cfreq_prof(1:7));
        [fitted_curve_l, goodness] = fit(-xax(2:end), contact_freq(2:end), ft, opts);
        plot(xax_o/30*(13-7)+7,fitted_curve_r(xax_o),'b-')
        plot(-xax_o/30*(13-7)+7,fitted_curve_l(xax_o),'b-')
        disp('dox');
        fprintf('Decay constant: %.3f\n', fitted_curve_r.b)
        fprintf('Decay constant: %.3f\n', fitted_curve_l.b)
        plawfit_b(ii*2,1) = fitted_curve_l.b;
        plawfit_b(ii*2,2) = fitted_curve_r.b;
        
    ylim([0.1,0.59]);
%     ylim([0,1]);
    xlim([0.5,13.5]);
    set(gca,'XTick',[1,4,7,10,13],'XTickLabel',[-30,-15,0,15,30],'XColor','k','YColor','k');
%     set(gca,'XTick',[1,4,7,10,13],'XColor','k','YColor','k');
    xlabel('kb');
    ylabel('contact frequency');
    title(dataselected{vizidx(ii)});
end
% saveas(fh,[figDIR,'cfreq_memoryseries_plaw_fit.pdf']);
% saveas(fh,[figDIR,'cfreq_memoryseries_ctrl.pdf']);

%% bar plot, median Rg
% figDIR = '/Users/tfuji/Documents/tempanalysis/_figures_v3/fig4/';

cmap = [0.3,0.3,0.8; 0.7,0.7,0.9];
idxK = [1, 2];
idxH = [9, 10];
idxD = [12, 13];
idx_st = 1; idx_ed = 13;
close all;
for ii = 8%1:7
    fh = SizedFig(15,25);
    hold on;
    switch ii
        case 1
            idx = idxK; 
            cond_names = 'KRAB-1day';
            xticklab = {'day0','day1'};
        case 2
            idx = idxH; 
            cond_names = 'HDAC4-1day';
            xticklab = {'day0','day1'};
        case 3
            idx = idxD; % idxKm; 
            cond_names = 'DNMT3B-irr'; % 'KRAB147';
            xticklab = {'no dox','irreversible'};
        case 4
            idx = [9, 11]; % idxKm; 
            cond_names = 'HDAC4-5days'; % 'KRAB147';
            xticklab = {'day0','day5'};
        case 5
            idx = [1,4]; 
            cond_names = 'KRAB-irr';
            xticklab = {'no dox','irreversible'};
        case 6
            idx = [1,5]; 
            cond_names = 'KRAB-react';
            xticklab = {'no dox','reactivated'};
        case 7
            idx = [1,3]; 
            cond_names = 'KRAB-5days';
            xticklab = {'day0','day5'};
        case 8
            idx = [17,18]; 
            cond_names = 'VP64-1day';
            xticklab = {'day0','day1'};
    end

    d0 = reshape(medmats(idx(1)).Rg(idx_st,idx_ed,:),[],1);
    d1 = reshape(medmats(idx(2)).Rg(idx_st,idx_ed,:),[],1);
%     d5 = reshape(medmats(idx(3)).Rg(idx_st,idx_ed,:),[],1);
%     Rg_series = [d0,d1,d5];
    Rg_tmp = [d0,d1];
    valid_idx = sum(~isnan(Rg_tmp),2) == 2;
    Rg_tmp = Rg_tmp(valid_idx,:);
    Rg_series_avg = mean(Rg_tmp,1);

    bar(Rg_series_avg,'FaceColor',cmap(1,:));
    plot(Rg_tmp','o-','MarkerFaceColor','w','Color',...
        cmap(2,:),'LineWidth',1);
    ylim([152,205]);
    ylim([152,220]); % for VP64
    xlim([0.2,2.8]);
%     set(gca,'XTick',1:3,'XTickLabel',{'day0','day1','day5'},...
%         'XColor','k','YColor','k');
    set(gca,'XTick',1:2,'XTickLabel',xticklab,...
        'XColor','k','YColor','k');
    ylabel('median Rg (-30:30kb) (nm)');
    [~,pt1] = ttest(Rg_tmp(:,1),Rg_tmp(:,2));
%     [~,pt5] = ttest(Rg_series(:,1),Rg_series(:,3));
    disp([pt1]);
%     saveas(fh,[figDIR,'bar_medianRg_',cond_names,'.pdf']);
end


%%
writematrix(Rg_series,[figDIR,'Rg_series.csv']);

%%
%% Mod vs Rg correlation
%%
%%
cnrDIR = '/Users/tfuji/Documents/tempanalysis/CUTnRUN/CNR_chr19/';
load([cnrDIR,'H3K9me3_KHDKmut_1kbbin.mat']);
st_kb = -30;
ed_kb = 30;
% st_kb = -40;
% ed_kb = 40;
pltrange = find(xax == st_kb):find(xax == ed_kb);
markprof_selected_avg_me = mean(markprof_selected(pltrange,:),1);
markprof_selected_avg_me = [markprof_selected_avg_me(1:2:end);markprof_selected_avg_me(2:2:end)];
% fill NaN for HDAC4-1day since I didn't measure
markprof_selected_avg_me = [markprof_selected_avg_me(:,1:9),nan(size(markprof_selected_avg_me,1),1),markprof_selected_avg_me(:,10:12)];
%
load([cnrDIR,'H3K9ac_KHDKmut_100bpbin.mat']);
st_kb = -10;
ed_kb = 10;
pltrange = find(xax == st_kb):find(xax == ed_kb);
pltrange = pltrange(pltrange < pEFregion(1) | pltrange > pEFregion(end));
markprof_selected_avg_ac = mean(markprof_selected(pltrange,:),1);
markprof_selected_avg_ac = [markprof_selected_avg_ac(1:2:end);markprof_selected_avg_ac(2:2:end)];


%% Mod vs med Rg
% close all
% fitfunc = @(b,x)(b(1)*log(x) + b(2));
fitfunc = @(b,x)(b(1)*x + b(2));
% AMratio = mean(markprof_selected_avg_ac,1)./mean(markprof_selected_avg_me,1);
Ac_avg = mean(markprof_selected_avg_ac,1);
Me_avg = mean(markprof_selected_avg_me,1);
axislim = {[0, 14, 172,202],[-0.5, 8, 172,202]};
% axislim = {[0, 14, 0.9,1.05],[-0.5, 10, 0.9,1.05]};
xlabs = {'H3K9ac','H3K9me3'};
Rg_median_avg = nanmean(Rg_series(1:2,:),1);
Rg_median_sem = nanstd(Rg_series(1:2,:),1)./sqrt(sum(~isnan(Rg_series(1:2,:)),1));
ydata = Rg_median_avg;
ysem = Rg_median_sem;
markershape = {'o','o','square','x','o'};
markertype = {'bo','bo','bo','bo','bo'};
markertypeerr = {'b-','b-','b-','b-','b-'};
conditiongrouped = {[1,2,3],[6,7,8],[9,10,11],[12,13],[4,5]};
% conditiongrouped = {[1,2,3],[6,7,8],[10],[13],[4],[5]};
for ii = 2
    switch ii
        case 1
            xdata = Ac_avg;
            xerr = max(markprof_selected_avg_ac,[],1)-Ac_avg;
        case 2
            xdata = Me_avg;
            xerr = max(markprof_selected_avg_me,[],1)-Me_avg;
    end
    
    fh = SizedFig(16,30);
    hold on;
%     plot(xax,yax,'k-');
    for jj = 1:length(markershape)
        idx = conditiongrouped{jj};
        plot(xdata(idx), ydata(idx),markertype{jj},'MarkerSize',15,'LineWidth',2,'Marker',markershape{jj});
        plot([xdata(idx);xdata(idx)],[ydata(idx)-ysem(idx);ydata(idx)+ysem(idx)],markertypeerr{jj});
        plot([xdata(idx)-xerr(idx);xdata(idx)+xerr(idx)],[ydata(idx);ydata(idx)],markertypeerr{jj});
    end
    axis(axislim{ii});
    ylabel('median Rg (nm)');
    xlabel(xlabs{ii});
%     set(gca,'XScale','log');
    set(gca,'FontSize',15,'XColor','k','YColor','k');   
    if ii == 2
%         saveas(fh,[figDIR,'Rg_vs_H3K9me3_memoryseries.pdf'])
%         saveas(fh,[figDIR,'Rg_vs_H3K9me3_nonKRABrec.pdf'])
    end
end
%% import RNAFISH, memory series
parDIR = '/Users/tfuji/Documents/tempanalysis/_figures_v3/fig4_5/A6RNAFISH_data/';
davgrna = [];
davgrna.spotnum = [];
davgrna.silfrac = [];
davgrna.dname = {};
for kk = 1:2
    switch kk
       case 1
           load([parDIR,'A6-RNAFISH_AvgObjByCondition_20220802.mat']);
           df1 = AvgObjsByCondition;
       case 2
           load([parDIR,'A6-RNAFISH_AvgObjByCondition_20220805.mat']);
           df1 = AvgObjsByCondition;
    end
    drange = 1:length(df1);
    for ii = 1:length(drange)

        criteria = (df1(ii).area > 7500);
        tmp = df1(ii).cy5spots(criteria);
        davgrna.spotnum = [davgrna.spotnum;mean(tmp)];
        davgrna.silfrac = [davgrna.silfrac;sum(tmp < 30)/length(tmp)];
        davgrna.dname = [davgrna.dname,df1(ii).dname];

    end
end
% mrna = (davgrna.silfrac([1,4,5,6,2,3,10,11,12])+davgrna.silfrac([1,4,5,6,2,3,10,11,12]+12))/2;
mrna = (davgrna.silfrac([1,4,5,6,2,3])+davgrna.silfrac([1,4,5,6,2,3]+12))/2;

npyDIR = '/Users/tfuji/Library/CloudStorage/GoogleDrive-tfuji@stanford.edu/My Drive/flow_analysis/rawdata_KRAB/KRAB_A6vsB9_OFFfracDay16_tgt-dox-rec-rep.npy';
test = readNPY(npyDIR);
memfrac = [reshape(test(4,1,2,:),[],1), reshape(test(5,1,1,:),[],1),...
    reshape(test(5,2,1,:),[],1), reshape(test(5,2,2,:),[],1),...
    reshape(test(4,2,1,:),[],1),reshape(test(4,2,2,:),[],1)];
memfrac = sqrt(memfrac);
%% gene expression vs med Rg
fitfunc = @(b,x)(b(1)*x + b(2));
axislim = {[-0.1, 1.1, 172,202],[-0.1, 0.8, 175,201]};
% axislim = {[-0.1, 1.1, 172,202],[-0.1, 0.8, 158,202]};
% axislim = {[-0.1, 1.1, 172,202],[-0.1, 0.8, 0.84,1.02]};
xlabs = {'Slienced fraction','memory fraction'};
markershape = {'o','o','o'};
markertype = {'bo','bo','bo'};
markertypeerr = {'b-','b-','b-'};
conditiongrouped = {[1,5,6],[3,4]};
% conditiongrouped = {[1,5,6],[2,3,4],[7,8,9]};
% conditiongrouped = {[1,5,6],[2,3,4],[10,11,12]};
for ii = 2
    switch ii
        case 1
            xdata = mrna';
            ydata = Rg_series(:,[1,6,7,8,2,3]);
%             ydata = Me_avg([1,6,7,8,2,3])';
        case 2
            xdata = mean(memfrac,1);
            xerr = std(memfrac,1);%/sqrt(size(memfrac,1));
            ydata = Rg_series(1:2,[1,6,7,8,2,3]);
%             ydata = ydata./repmat(ydata(:,1),[1 size(ydata,2)]);
    end
    ydata_avg = nanmean(ydata,1);
    ydata_sem = nanstd(ydata,1);%./sqrt(sum(~isnan(ydata),1));
    valid_indices = ~isnan(xdata) & ~isnan(ydata_avg);
    xdata = xdata(valid_indices);
    ydata_avg = ydata_avg(valid_indices);
    beta0 = [2;180];
    beta = nlinfit(xdata,ydata_avg,fitfunc,beta0);
    xax = linspace(axislim{ii}(1),axislim{ii}(2),100);
    yax = fitfunc(beta,xax);
    corrcoef(xdata, ydata_avg)

    fh = SizedFig(16,30);
    hold on;
    plot(xax,yax,'k-');
    for jj = 1:length(conditiongrouped)
        idx = conditiongrouped{jj};
        plot(xdata(idx), ydata_avg(idx),markertype{jj},'MarkerSize',15,'LineWidth',2,'Marker',markershape{jj});
        plot([xdata(idx);xdata(idx)],[ydata_avg(idx)-ydata_sem(idx);ydata_avg(idx)+ydata_sem(idx)],markertypeerr{jj},'LineWidth',1);
        plot([xdata(idx)-xerr(idx);xdata(idx)+xerr(idx)],[ydata_avg(idx);ydata_avg(idx)],markertypeerr{jj},'LineWidth',1);
%         for kk = 1:length(idx)
%             plot(xdata(:,idx(kk)),ydata(:,idx(kk)),'r.','MarkerSize',5);
%         end
    end
    axis(axislim{ii});
    ylabel('median Rg (nm)');
%     ylabel('median Rg ratio (KRAB nodox=1)');
    xlabel(xlabs{ii});
    set(gca,'FontSize',15,'XColor','k','YColor','k');   
    switch ii
        case 2
%         saveas(fh,[figDIR,'RatiomedRg_vs_memory_allrep_1+4points_willeachrepdots.pdf']);
%         saveas(fh,[figDIR,'Rg_vs_memory_allrep_1+4points_willeachrepdots.pdf']);
%         saveas(fh,[figDIR,'Rg_vs_memory_allrep_1+4points.pdf']);
%         saveas(fh,[figDIR,'Rg_vs_memory_allrep_krabhdac.pdf']);
%         saveas(fh,[figDIR,'Rg_vs_memory_onlyrep12.pdf']);
    end
end
%% gene expression vs chrom mod
% fitfunc = @(b,x)(b(1)*log(x) + b(2));
fitfunc = @(b,x)(b(1)*x + b(2));
% AMratio = mean(markprof_selected_avg_ac,1)./mean(markprof_selected_avg_me,1);
Ac_avg = mean(markprof_selected_avg_ac,1);
Me_avg = mean(markprof_selected_avg_me,1);
axislim = {[-0.1, 0.8, 0,15],[-0.1, 0.8, 0,8]};
% axislim = {[-0.1, 1.1, 0,8],[-0.1, 0.8, 0,8]};
xlabs = {'memory fraction','memory fraction'};
ylabs = {'H3K9ac','H3K9me3'};
markershape = {'o','o','o'};
markertype = {'bo','bo','bo'};
markertypeerr = {'b-','b-','b-'};
conditiongrouped = {[1,5,6],[3,4]};
for ii = 1:2
    switch ii
        case 1
            xdata = mean(memfrac,1);
            xerr = std(memfrac,1);%/sqrt(size(memfrac,1));
            ydata = Ac_avg([1,6,7,8,2,3]);
            yerr = max(markprof_selected_avg_ac(:,[1,6,7,8,2,3]),[],1)-ydata;
        case 2
            xdata = mean(memfrac,1);
            xerr = std(memfrac,1);%/sqrt(size(memfrac,1));
            ydata = Me_avg([1,6,7,8,2,3]);
            yerr = max(markprof_selected_avg_me(:,[1,6,7,8,2,3]),[],1)-ydata;
    end
    ydata_avg = nanmean(ydata,1);
    ydata_sem = nanstd(ydata,1);%./sqrt(sum(~isnan(ydata),1));
    valid_indices = ~isnan(xdata) & ~isnan(ydata_avg);
    xdata = xdata(valid_indices);
    ydata_avg = ydata_avg(valid_indices);
    beta0 = [2;180];
    beta = nlinfit(xdata,ydata_avg,fitfunc,beta0);
    xax = linspace(axislim{ii}(1),axislim{ii}(2),100);
    yax = fitfunc(beta,xax);
    corrcoef(xdata, ydata_avg)

    fh = SizedFig(16,30);
    hold on;
    plot(xax,yax,'k-');
    for jj = 1:length(conditiongrouped)
        idx = conditiongrouped{jj};
        plot(xdata(idx), ydata_avg(idx),markertype{jj},'MarkerSize',15,'LineWidth',2,'Marker',markershape{jj});
        plot([xdata(idx);xdata(idx)],[ydata_avg(idx)-yerr(idx);ydata_avg(idx)+yerr(idx)],markertypeerr{jj},'LineWidth',1);
        plot([xdata(idx)-xerr(idx);xdata(idx)+xerr(idx)],[ydata_avg(idx);ydata_avg(idx)],markertypeerr{jj},'LineWidth',1);
    end
%     axis(axislim{ii});
    ylabel(ylabs{ii});
    xlabel(xlabs{ii});
    set(gca,'FontSize',15,'XColor','k','YColor','k');   
    switch ii
%         case 1
%         saveas(fh,[figDIR,'H3K9ac_vs_memory_allrep_1+4points.pdf']);
%         case 2
%         saveas(fh,[figDIR,'H3K9me3_vs_memory_allrep_1+4points.pdf']);
    end
end

%% gene expression vs med Rg; different marker for each condition
fitfunc = @(b,x)(b(1)*x + b(2));
axislim = {[-0.1, 1.1, 172,202],[-0.1, 0.8, 175,201]};
% axislim = {[-0.1, 1.1, 172,202],[-0.1, 0.8, 158,202]};
% axislim = {[-0.1, 1.1, 172,202],[-0.1, 0.8, 0.84,1.02]};
xlabs = {'Slienced fraction','memory fraction'};
markershape = {'o','s','o'};
markertype = {'bo','bs','bo'};
markertypeerr = {'b-','b-','b-'};
conditiongrouped = {[1,5,6],[2,3,4]};
% conditiongrouped = {[1,5,6],[2,3,4],[7,8,9]};
% conditiongrouped = {[1,5,6],[2,3,4],[10,11,12]};
color_eachmark = {'#0000FF','#0000FF'};
for ii = 2
    switch ii
        case 1
            xdata = mrna';
            ydata = Rg_series(:,[1,6,7,8,2,3]);
%             ydata = Me_avg([1,6,7,8,2,3])';
        case 2
            xdata = mean(memfrac,1);
            xerr = std(memfrac,1);%/sqrt(size(memfrac,1));
            ydata = Rg_series(1:2,[1,6,7,8,2,3]);
%             ydata = ydata./repmat(ydata(:,1),[1 size(ydata,2)]);
    end
    ydata_avg = nanmean(ydata,1);
    ydata_sem = nanstd(ydata,1);%./sqrt(sum(~isnan(ydata),1));
    valid_indices = ~isnan(xdata) & ~isnan(ydata_avg);
    xdata = xdata(valid_indices);
    ydata_avg = ydata_avg(valid_indices);
    beta0 = [2;180];
    mdl = fitnlm(xdata, ydata_avg,fitfunc,beta0);
    xax = linspace(axislim{ii}(1),axislim{ii}(2),100);
    yax = fitfunc(mdl.Coefficients.Estimate,xax);
    disp(mdl.Rsquared.Ordinary);
        
    fh = SizedFig(16,30);
    plot(xax,yax,'k-');
    hold on;
    for jj = 1:length(conditiongrouped)
        idx = conditiongrouped{jj};
        plot(xdata(idx), ydata_avg(idx),'--','MarkerSize',15,'LineWidth',1,'Marker',markershape{jj},'Color',hex2rgb(color_eachmark{ii})*jj/2);
        plot(xdata(idx), ydata_avg(idx),markertype{jj},'MarkerSize',15,'LineWidth',2,'Marker',markershape{jj},'Color',hex2rgb(color_eachmark{ii})*jj/2);
%         plot(xdata(idx), ydata_avg(idx),[markertype{jj},'-'],'MarkerSize',15,'LineWidth',2,'Marker',markershape{jj});
        plot([xdata(idx);xdata(idx)],[ydata_avg(idx)-ydata_sem(idx);ydata_avg(idx)+ydata_sem(idx)],markertypeerr{jj},'LineWidth',1,'Color',hex2rgb(color_eachmark{ii})*jj/2);
        plot([xdata(idx)-xerr(idx);xdata(idx)+xerr(idx)],[ydata_avg(idx);ydata_avg(idx)],markertypeerr{jj},'LineWidth',1,'Color',hex2rgb(color_eachmark{ii})*jj/2);
%         for kk = 1:length(idx)
%             plot(xdata(:,idx(kk)),ydata(:,idx(kk)),'r.','MarkerSize',5);
%         end
%         mdl = fitnlm(xdata(idx), ydata_avg(idx),fitfunc,beta0);
%         xax = linspace(axislim{ii}(1),axislim{ii}(2),100);
%         yax = fitfunc(mdl.Coefficients.Estimate,xax);
%         disp(mdl.Rsquared.Ordinary);
%         plot(xax,yax,'k-');
    end
    axis(axislim{ii});
    ylabel('median Rg (nm)');
%     ylabel('median Rg ratio (KRAB nodox=1)');
    xlabel(xlabs{ii});
    set(gca,'FontSize',15,'XColor','k','YColor','k');   
    box off;
    switch ii
        case 2
%         saveas(fh,[figDIR,'RatiomedRg_vs_memory_allrep_1+4points_willeachrepdots.pdf']);
%         saveas(fh,[figDIR,'Rg_vs_memory_allrep_1+4points_willeachrepdots.pdf']);
%         saveas(fh,[figDIR,'Rg_vs_memory_allrep_1+4points.pdf']);
%         saveas(fh,[figDIR,'Rg_vs_memory_allrep_2+4points.pdf']);
%         saveas(fh,[figDIR,'Rg_vs_memory_allrep_2+4points_wLine.pdf']);
%         saveas(fh,[figDIR,'Rg_vs_memory_allrep_krabhdac.pdf']);
%         saveas(fh,[figDIR,'Rg_vs_memory_onlyrep12.pdf']);
    end
end
%% gene expression vs med Rg; different marker for each condition, KRAB150
fitfunc = @(b,x)(b(1)*x + b(2));
axislim = {[-0.1, 1.1, 172,202],[-0.1, 0.8, 175,201]};
xlabs = {'Slienced fraction','memory fraction'};
markershape = {'^','s','o'};
markertype = {'b^','bs','bo'};
markertypeerr = {'b-','b-','b-'};
conditiongrouped = {[7,8,9]};
color_eachmark = {'#0000FF','#0000FF'};
for ii = 2
    switch ii
        case 1
            xdata = mrna';
            ydata = Rg_series(:,[1,6,7,8,2,3]);
%             ydata = Me_avg([1,6,7,8,2,3])';
        case 2
            xdata = [0,0.3,0.6,0,0.3,0.6,0,0.3,0.6];
            xerr = [0,0,0];
            ydata = Rg_series(1:2,[1,6,7,8,2,3,14,15,16]);
%             ydata = ydata./repmat(ydata(:,1),[1 size(ydata,2)]);
    end

        
    fh = SizedFig(16,30);
    hold on;
    for jj = 1:length(conditiongrouped)
        idx = conditiongrouped{jj};
        plot(xdata(idx), ydata_avg(idx),markertype{jj},'MarkerSize',15,'LineWidth',2,'Marker',markershape{jj},'Color',hex2rgb(color_eachmark{ii}));
        plot([xdata(idx);xdata(idx)],[ydata_avg(idx)-ydata_sem(idx);ydata_avg(idx)+ydata_sem(idx)],markertypeerr{jj},'LineWidth',1,'Color',hex2rgb(color_eachmark{ii}));
    end
    axis(axislim{ii});
    ylabel('median Rg (nm)');
    xlabel(xlabs{ii});
    set(gca,'FontSize',15,'XColor','k','YColor','k');   
    box off;
    switch ii
        case 2
%         saveas(fh,[figDIR,'Rg_vs_memory_allrep_KRAB150.pdf']);
    end
end
%% gene expression vs chrom mod; different marker for each condition
% fitfunc = @(b,x)(b(1)*log(x) + b(2));
fitfunc = @(b,x)(b(1)*x + b(2));
% AMratio = mean(markprof_selected_avg_ac,1)./mean(markprof_selected_avg_me,1);
Ac_avg = mean(markprof_selected_avg_ac,1);
Me_avg = mean(markprof_selected_avg_me,1);
axislim = {[-0.1, 0.8, 0,15],[-0.1, 0.8, 0,8]};
% axislim = {[-0.1, 1.1, 0,8],[-0.1, 0.8, 0,8]};
xlabs = {'memory fraction','memory fraction'};
ylabs = {'H3K9ac','H3K9me3'};
markershape = {'o','s','o'};
markertype = {'bo','bs','bo'};
markertypeerr = {'b-','b-','b-'};
conditiongrouped = {[1,5,6],[2,3,4]};
color_eachmark = {'#1AE61A','#CC66CC'};
for ii = 1:2
    switch ii
        case 1
            xdata = mean(memfrac,1);
            xerr = std(memfrac,1);%/sqrt(size(memfrac,1));
            ydata = Ac_avg([1,6,7,8,2,3]);
            yerr = max(markprof_selected_avg_ac(:,[1,6,7,8,2,3]),[],1)-ydata;
        case 2
            xdata = mean(memfrac,1);
            xerr = std(memfrac,1);%/sqrt(size(memfrac,1));
            ydata = Me_avg([1,6,7,8,2,3]);
            yerr = max(markprof_selected_avg_me(:,[1,6,7,8,2,3]),[],1)-ydata;
    end
    ydata_avg = nanmean(ydata,1);
    ydata_sem = nanstd(ydata,1);%./sqrt(sum(~isnan(ydata),1));
    valid_indices = ~isnan(xdata) & ~isnan(ydata_avg);
    xdata = xdata(valid_indices);
    ydata_avg = ydata_avg(valid_indices);
    beta0 = [2;180];
    mdl = fitnlm(xdata, ydata_avg,fitfunc,beta0);
    xax = linspace(axislim{ii}(1),axislim{ii}(2),100);
    yax = fitfunc(mdl.Coefficients.Estimate,xax);
    disp(mdl.Rsquared.Ordinary);
    
    fh = SizedFig(16,30);
    hold on;
    plot(xax,yax,'k-');
    for jj = 1:length(conditiongrouped)
        idx = conditiongrouped{jj};
        plot(xdata(idx), ydata_avg(idx),'--','MarkerSize',15,'LineWidth',1,'Marker',markershape{jj},'Color',hex2rgb(color_eachmark{ii})*jj/2);
        plot(xdata(idx), ydata_avg(idx),markertype{jj},'MarkerSize',15,'LineWidth',2,'Marker',markershape{jj},'Color',hex2rgb(color_eachmark{ii})*jj/2);
        plot([xdata(idx);xdata(idx)],[ydata_avg(idx)-yerr(idx);ydata_avg(idx)+yerr(idx)],markertypeerr{jj},'LineWidth',1,'Color',hex2rgb(color_eachmark{ii})*jj/2);
        plot([xdata(idx)-xerr(idx);xdata(idx)+xerr(idx)],[ydata_avg(idx);ydata_avg(idx)],markertypeerr{jj},'LineWidth',1,'Color',hex2rgb(color_eachmark{ii})*jj/2);
        
%         mdl = fitnlm(xdata(idx), ydata_avg(idx),fitfunc,beta0);
%         xax = linspace(axislim{ii}(1),axislim{ii}(2),100);
%         yax = fitfunc(mdl.Coefficients.Estimate,xax);
%         disp(mdl.Rsquared.Ordinary);
%         plot(xax,yax,'k-');
    end
    axis(axislim{ii});
    ylabel(ylabs{ii});
    xlabel(xlabs{ii});
    set(gca,'FontSize',15,'XColor','k','YColor','k');   
    switch ii
%         case 1
%         saveas(fh,[figDIR,'H3K9ac_vs_memory_allrep_2+4points_wLine.pdf']);
%         case 2
%         saveas(fh,[figDIR,'H3K9me3_vs_memory_allrep_2+4points_wLine.pdf']);
    end
end
%% correlation coefficient matrix
{'Silencing','Memory','K9ac','K9me3','Rg'};
Rg_median_avg = nanmean(Rg_series(1:2,:),1);
% Tout_Mod_Rg = [Ac_avg([1,6,7,8,2,3])',Me_avg([1,6,7,8,2,3])',Rg_median_avg([1,6,7,8,2,3])',mrna,mean(memfrac,1)'];
Tout_Mod_Rg = [Ac_avg([1,7,8,2,3])',Me_avg([1,7,8,2,3])',Rg_median_avg([1,7,8,2,3])',mrna([1,3,4,5,6]),mean(memfrac(:,[1,3,4,5,6]),1)'];
cf = nan(size(Tout_Mod_Rg,2),size(Tout_Mod_Rg,2));
for ii = 1:(size(cf,1))
    for jj = (ii):size(cf,1)
%         tmp = corrcoef(Tout_Mod_Rg(:,ii),Tout_Mod_Rg(:,jj));
%         cf(ii,jj) = tmp(1,2);
        beta0 = [2;180];
        mdl = fitnlm(Tout_Mod_Rg(:,ii),Tout_Mod_Rg(:,jj),fitfunc,beta0);
        cf(ii,jj) = mdl.Rsquared.Ordinary;
%         mdl = corrcoef(Tout_Mod_Rg(2:end,ii),Tout_Mod_Rg(2:end,jj));
%         cf(ii,jj) = mdl(1,2);
    end
end
fh = SizedFig(25,30);
imagesc(abs(cf(1:3,4:5))');
axis image;
caxis([0.7,0.8]);
% caxis([0.7,0.9]);
set(gca,'XTick',[],'YTick',[]);
box off;
colorbar;
% saveas(fh,[figDIR,'corrcoef_Y=sil-mem_X=Ac-Me-Rg.pdf']);
% saveas(fh,[figDIR,'corrcoef_R2_Y=sil-mem_X=Ac-Me-Rg.pdf']);
%% bar memory series; Rg
figDIR = '/Users/tfuji/Documents/tempanalysis/_figures_v3/fig4/';

cmap = [0.3,0.3,0.8; 0.7,0.7,0.9];
idx = [1,7, 8, 2, 3];
idx_st = 1; idx_ed = 13;
close all;
fh = SizedFig(15,18);
%     subplot(1,2,ii);
hold on;
% Rg_series_tmp = zeros(size(medmats(1).Rg,3),length(idx));
Rg_series_tmp = zeros(2,length(idx));
for ii = 1:length(idx)
    Rg_series_tmp(:,ii) = reshape(medmats(idx(ii)).Rg(idx_st,idx_ed,1:2),[],1);
end
Rg_series_avg = nanmean(Rg_series_tmp,1);
Rg_series_sem = nanstd(Rg_series_tmp,1)./sqrt(sum(~isnan(Rg_series_tmp),1));

bar(Rg_series_avg,'FaceColor',cmap(1,:));
% plot([1:length(idx);1:length(idx)],[Rg_series_avg-Rg_series_sem;Rg_series_avg+Rg_series_sem],'k-','LineWidth',1);
plot([1:length(idx);1:length(idx)],Rg_series_tmp,'ko','LineWidth',1,'MarkerFaceColor','w');
% plot(Rg_series','o','MarkerFaceColor','w','Color',...
%     cmap(2,:),'LineWidth',1);
ylim([172,201]);
xlim([0.3,5.7]);
set(gca,'XColor','k','YColor','k');
ylabel('median Rg (-30:30kb) (nm)');
% saveas(fh,[figDIR,'bar_memseries_Rg.pdf']);

%% bar memory series; H3K9ac
figDIR = '/Users/tfuji/Documents/tempanalysis/_figures_v3/fig4/';
idx = [1,7, 8, 2, 3];
Ac_avg = mean(markprof_selected_avg_ac,1);
close all;
fh = SizedFig(15,18);
%     subplot(1,2,ii);
hold on;
ydata = Ac_avg(idx);
cmap = [0.1, 0.9, 0.1];
bar(ydata,'FaceColor',cmap);
plot([1:length(idx);1:length(idx)],markprof_selected_avg_ac(:,idx),'ko','LineWidth',1,'MarkerFaceColor','w');
% ylim([175,200]);
xlim([0.3,5.7]);
set(gca,'XColor','k','YColor','k');
ylabel('H3K9ac(-10:10kb) (reads/kb)');
% saveas(fh,[figDIR,'bar_memseries_H3K9ac.pdf']);
%% bar memory series, H3K9me3
figDIR = '/Users/tfuji/Documents/tempanalysis/_figures_v3/fig4/';
idx = [1,7, 8, 2, 3];
close all;
Me_avg = mean(markprof_selected_avg_me,1);
fh = SizedFig(15,18);
%     subplot(1,2,ii);
hold on;
ydata = Me_avg(idx);
cmap = [0.8, 0.4, 0.8]; % ; 0.9, 0.5, 0.9; 0.9, 0.7, 0.9];
bar(ydata,'FaceColor',cmap);
plot([1:length(idx);1:length(idx)],markprof_selected_avg_me(:,idx),'ko','LineWidth',1,'MarkerFaceColor','w');
xlim([0.3,5.7]);
set(gca,'XColor','k','YColor','k');
ylabel('H3K9me3(-20:40kb) (reads/kb)');
% saveas(fh,[figDIR,'bar_memseries_H3K9me3.pdf']);

%%

%%
%%
%%
%% combined all reps to calculate Rg
%%
%%
idx_st = 1; idx_ed = 13;
% idx_st = 5; idx_ed = 9;
Rg_combined = {};
for ii = 1:length(dataselected)
    tmpRg = reshape(dfig2.Rg(idx_st,idx_ed,strcmp(dfig2.dname,dataselected(ii))),[],1);
    Rg_combined = [Rg_combined,tmpRg];
end


%% mod vs Rg
AMratio = mean(markprof_selected_avg_ac,1)./mean(markprof_selected_avg_me,1);
Ac_avg = mean(markprof_selected_avg_ac,1);
Me_avg = mean(markprof_selected_avg_me,1);
Rg_median = cellfun(@median,Rg_combined);
fitfunc = @(b,x)(b(1)*log(x) + b(2));
% axislim = {[-1, 20, 180,200],[-0.5, 4.7, 180,200],[0.3, 299, 180,200]};
axislim = {[0.5, 14, 180,200],[0.01, 10, 180,200],[0.1, exp(6), 180,200]};
xlabs = {'H3K9ac','H3K9me3','H3K9ac/H3K9me3'};
SizedFig(40,25);
for ii = 1:3
    subplot(1,3,ii);
    hold on;
    switch ii
        case 1
            xdat = Ac_avg;
            ydat = Rg_median;
        case 2
            xdat = Me_avg;
            ydat = Rg_median;
        case 3
            xdat = Ac_avg./Me_avg;
            ydat = Rg_median;
    end
    valid_idx = ~isnan(xdat) & ~isnan(ydat);
    xdat = xdat(valid_idx(1:end)); ydat = ydat(valid_idx(1:end));
    beta0 = [2;180];
    beta = nlinfit(xdat,ydat,fitfunc,beta0);
%     xax = linspace(0,axislim{ii}(2),2000);
    xax = linspace(axislim{ii}(1),axislim{ii}(2),2000);
    yax = fitfunc(beta,xax);
    plot(xax,yax,'k-');
    plot(xdat,ydat,'bo','MarkerSize',8,'LineWidth',2,...
        'MarkerFaceColor','w');
    axis(axislim{ii});
    if ii == 3 % || ii == 2 || ii == 1
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    set(gca,'XTick',[10^0,10^1,10^2]);
    end
    ylabel('median Rg (nm)');
    xlabel(xlabs{ii});
    set(gca,'FontSize',10);
end



%% recalc Rg

for objidx = 1:length(Objs)
    df = Objs(objidx).df;

    % --- filter >50% positive rate
    for i = 1:length(df)
        tmp = df(i).coordfilt;
        q1 = reshape(sum(~isnan(tmp(1:19,1,:)),1),[],1);
        posidx = (q1 >= 10);
        df(i).coordfilt50p = df(i).coordfilt(:,:,posidx);
        df(i).dmatfilt50p = df(i).dmatfilt(:,:,posidx);
    end

    % --- linear interpolated matrix
    for i  = 1:length(df)
        tmpd = df(i).coordfilt50p(1:end-1,:,:);
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
        nanidx = reshape(isnan(sum(sum(dmatinterp,1),2)),[],1);
        df(i).dmatinterp = dmatinterp(:,:,~nanidx);
        df(i).coordinterp = coordinterp(:,:,~nanidx);
    end

    % --- pick up upper triangle for dimentionality reduction
    idxrange = 4:16;
%     idxrange = 1:13;
    for ii = 1:length(df)
        allpairs = [];
        tmp = df(ii).dmatinterp(idxrange,idxrange,:);
        tmp_triu = ~triu(ones(size(tmp,1),size(tmp,2)))';
        for i = 1:size(tmp,3)
            tmpmat = tmp(:,:,i);
            allpairs = cat(1,allpairs,tmpmat(tmp_triu)');
        end
        df(ii).pairs = allpairs;
    end
    
    Objs(objidx).df = df;
end




%% -----------------------------------------------------------------
%% -------- median distance map for each bioreps and ranksum -------
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
        tmp = dfig2.dmat(:,:,strcmp(dfig2.dname,dataselected(ii)) & dfig2.rep == jj);
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
        tmp = dfig2.dmatfilt(:,:,strcmp(dfig2.dname_dmatfilt,dataselected(ii)) & dfig2.rep_dmatfilt == jj);
        tmp = nanmedian(tmp,3);
        medmat2 = cat(3,medmat2,tmp);
        tmp = dfig2.Rg(:,:,strcmp(dfig2.dname,dataselected(ii)) & dfig2.rep == jj);
        tmp = nanmedian(tmp,3);
        medmat3 = cat(3,medmat3,tmp);

        % pval
        tmp1 = dfig2.Rg(1,end,strcmp(dfig2.dname,dataselected(ii)) & dfig2.rep == jj);
        tmp2 = dfig2.Rg(1,end,strcmp(dfig2.dname,dataselected(1)) & dfig2.rep == jj);
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
% box on;