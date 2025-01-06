%%
parDIR = '/Users/tfuji/Documents/tempanalysis/EM-seq/';
%% params
DIR1 = '/Users/tfuji/Documents/tempanalysis/EM-seq/20230728_EM-seq_B9-KRAB_5d_irr_byAbby/heatmaps/';
fnames = dir(DIR1);
fnames = fnames(3:end);
dnames = {};
for ii = 1:length(fnames)
    tmpcell = split(fnames(ii).name,'.');
    if strcmp(tmpcell{end},'matrix') && strcmp(tmpcell{end-1},'clustered')
        dnames = [dnames,fnames(ii).name];
        disp([num2str(length(dnames)),':',fnames(ii).name]);
    end
end
%% params
DIR2 =  '/Users/tfuji/Documents/tempanalysis/EM-seq/20231013_EM-seq_B9-KRAB_5d_irr_react/heatmaps/';
fnames = dir(DIR2);
fnames = fnames(3:end);
dnames2 = {};
for ii = 1:length(fnames)
    tmpcell = split(fnames(ii).name,'.');
    if strcmp(tmpcell{end},'matrix') && strcmp(tmpcell{end-1},'clustered')
        dnames2 = [dnames2,fnames(ii).name];
        disp([num2str(length(dnames2)),':',fnames(ii).name]);
    end
end
%% params
DIR3 =  '/Users/tfuji/Documents/tempanalysis/EM-seq/20240317_EM-seq_HDAC4_timecourse_DNMTirr/Matrices_2/';
fnames = dir(DIR3);
fnames = fnames(4:end);
fnames = fnames(1:16);
dnames3 = {};
for ii = 1:length(fnames)
    tmpcell = split(fnames(ii).name,'.');
    if strcmp(tmpcell{end},'matrix') && strcmp(tmpcell{end-1},'full_unclustered')
        dnames3 = [dnames3,fnames(ii).name];
        disp([num2str(length(dnames3)),':',fnames(ii).name]);
    end
end
% L166	2-1
% L167	5-1
% L168	10-1
% L169	2-2
% L170	5-2
% L171	10-2
% L172	H+1
% L173	H+2
% L174	H-1
% L175	H-2
% L176	D-1
% L177	D-2
% L178	H01
% L179	H02
% L180	H51
% L181	H52

%%
dname = {'DNMT3B no dox','DNMT3B dox 12 days','DNMT3B/A6 no dox','DNMT3B/A6 dox 12 days',...
    'KRAB no dox rep1','KRAB dox 5 days rep1','irreversible rep1',...
    'KRAB no dox rep2','KRAB dox 5 days rep2','irreversible rep2',...
    'Abby control','reactivated rep1','reactivated rep2',...
    'KRAB dox5 release2 rep1','KRAB dox5 release5 rep1','KRAB dox5 release10 rep1',...
    'KRAB dox5 release2 rep2','KRAB dox5 release5 rep2','KRAB dox5 release10 rep2',...
    'HDAC4 reactivated rep1','HDAC4 reactivated rep2','HDAC4 irreversible rep1','HDAC4 irreversible rep2',...
    'DNMT3B irreversible rep1','DNMT3B irreversible rep2',...
    'HDAC4 no dox rep1','HDAC4 no dox rep2','HDAC4 dox 5days rep1','HDAC4 dox 5days rep2'};
dir1idx = [1,2,3,4,5,6,7,8,9, 10,11,NaN,NaN,nan(1,16)];
dir2idx = [1,2,3,4,5,6,7,8,9,NaN,10, 11, 12,nan(1,16)];
dir3idx = [nan(1,13),1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
%% import EM-seq results, heatmap met

% -- 1:unmet, 0:met -- %

cpgidx = [31  ,38  ,51  ,72  ,84  ,103 ,105 ,129 ,142 ,153 ,170 ,189 ,192 ,198 ,208 ,213 ,222 ,243 ,256 ,258 ,275 ,291 ,324 ,339 ,346 ,371 ,381 ,398 ,403 ,430 ,439 ,442 ,444 ,448 ,466 ,468 ,477 ,486 ,527 ,530 ,563 ,590 ,603 ,605 ,609 ,612 ,615 ,622 ,626 ,634 ,644 ,647 ,652 ,661 ,665 ,667 ,674 ,681 ,685 ,706 ,729 ,731 ,734 ,737 ,744 ,749 ,758]+1;

dnamet = {};
for ii = 1:length(dir1idx)
    tmp1 = [];
    if ~isnan(dir1idx(ii))
        tmp1 = readmatrix([DIR1,dnames{dir1idx(ii)}] ,'FileType', 'text');
        tmp1 = tmp1(2:end,2:end);
        tmp1 = tmp1(:,cpgidx([1:26,27,28,29,30,(end-26):end]));
        tmp1(:,[27,28,29,30]) = 2;
    end
    tmp2 = [];
    if ~isnan(dir2idx(ii))
        tmp2 = readmatrix([DIR2,dnames2{dir2idx(ii)}] ,'FileType', 'text');
        tmp2 = tmp2(2:end,2:end);
        tmp2 = tmp2(:,cpgidx([1:26,27,28,29,30,(end-26):end]));
        tmp2(:,[27,28,29,30]) = 2;
    end
    tmp3 = [];
    if ~isnan(dir3idx(ii))
        tmp3 = readmatrix([DIR3,dnames3{dir3idx(ii)}] ,'FileType', 'text');
        tmp3 = tmp3(2:end,2:end);
        tmp3 = tmp3(:,cpgidx([1:26,27,28,29,30,(end-26):end]));
        tmp3(:,[27,28,29,30]) = 2;
    end
    tmp = cat(1,tmp1,tmp2,tmp3);
    disp([size(tmp1,1),size(tmp2,1),size(tmp3,1),size(tmp,1)]);
    dnamet = [dnamet,tmp];
end

%% heatmap, xax adjusted
% actual index is like 50 -> 36 (-14 at the middle) -> 40 (+4 for spacing)
xaxlab = [10,20,50,60];
xax    = [10,20,40,50];
cmap = [1.0,1.0,1.0 ; 0.0, 0.0, 0.0 ; 0.9, 0.9, 0.9; 0.5,0,0];
SizedFig(15,100);
% visidx = 1:19;
visidx = [5,6,14,15,16,7,12]; % KRAB recruitment rep1
% visidx = [8,9,17,18,19,10,13]; % KRAB recruitment rep2
visidx = [26,27]; % KRAB recruitment rep1
for ii = 1:length(visidx)
    subplot(length(visidx),1,ii);
%     if ii <= 9
%     subplot(length(visidx),1,ii);
%     else
%     subplot(length(visidx),1,ii+1);
%     end
    tmp = dnamet{visidx(ii)};
    imagesc(tmp);
    colormap(cmap);
%     xlim([0,100]);
    caxis([-1,2]);
    title(dname{visidx(ii)});
    set(gca,'XTick',xax,'XTickLabel',xaxlab);
    set(gca,'XColor','k','YColor','k');
end

%% heatmap, xax adjusted; rep combined and sorted
% actual index is like 50 -> 36 (-14 at the middle) -> 40 (+4 for spacing)
xaxlab = [10,20,50,60];
xax    = [10,20,40,50];
cmap = [1.0,1.0,1.0 ; 0.0, 0.0, 0.0 ; 0.9, 0.9, 0.9; 0.5,0,0];
SizedFig(15,100);
% visidx = 1:19;
visidx = [5,6,14,15,16,7,12]; % KRAB recruitment rep1
visidx2 = [8,9,17,18,19,10,13]; % KRAB recruitment rep2
for ii = 1:length(visidx)
    subplot(length(visidx),1,ii);
    tmp = [dnamet{visidx(ii)};dnamet{visidx2(ii)}];
    [~,Idx] = sort(sum(tmp == 1,2)./(sum(tmp == 1,2)+sum(tmp == 0,2)),'descend');
    imagesc(tmp(Idx,:));
    colormap(cmap);
%     xlim([0,100]);
    caxis([-1,2]);
    title(dname{visidx(ii)});
    set(gca,'XTick',xax,'XTickLabel',xaxlab);
    set(gca,'XColor','k','YColor','k');
end


%% median of DNAmet fraction
MedDNAmet = [];
SizedFig(15,100);
visidx = [5,6,14,15,16,7,12];
visidx = [8,9,17,18,19,10,13];
visidx = 1:19;
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
MedDNAmet = cat(1,MedDNAmet([5,6,14,15,16,7,12]),MedDNAmet([8,9,17,18,19,10,13]));
DNAmetscore = nanmean(MedDNAmet,1);

% figure;
% plot(tmp,'bo');

%% plot average time course; including day 0
ymax = 0.75;
SizedFig(20,20);
hold on;
plot([5,5],[0,ymax],'k--','LineWidth',1,'Color',[0.9,0.2,0.2]);

plot([0,5,7,10,15],DNAmetscore(1:end-2),'r-','LineWidth',2,'Color',[0.9,0.2,0.2]);
plot([0,5,7,10,15],MedDNAmet(:,1:end-2),'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','w','Color',[0.9,0.2,0.2]);

plot([19,22],MedDNAmet(:,[6,7])','ro','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','w','Color',[0.9,0.2,0.2]);


text(19,max(MedDNAmet(:,6))+ymax/10,'irr.','HorizontalAlignment','center');
text(22,max(MedDNAmet(:,7))+ymax/10,'react.','HorizontalAlignment','center');
set(gca,'XTick',[0,5,7,10,15],'XColor','k','YColor','k');

title('DNA metylation');
xlabel('days');
ylabel('median DNAme fraction');
%% plot average time course; excluding day 0
ymax = 0.75;
fh = SizedFig(15,20);
hold on;
% plot([5,5],[0,ymax],'k--','LineWidth',1,'Color',[0.9,0.2,0.2]);

plot([15,19]-5,[DNAmetscore(end-2),mean(MedDNAmet(:,6))],'r--','LineWidth',2,'Color',[0.9,0.2,0.2]);

plot([5,7,10,15]-5,DNAmetscore(2:end-2),'r-','LineWidth',2,'Color',[0.9,0.2,0.2]);
plot([5,7,10,15]-5,MedDNAmet(:,2:end-2),'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','w','Color',[0.9,0.2,0.2]);

plot([19,22]-5,MedDNAmet(:,[6,7])','ro','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','w','Color',[0.9,0.2,0.2]);


text(19-5,max(MedDNAmet(:,6))+ymax/10,'irr.','HorizontalAlignment','center');
text(22-5,max(MedDNAmet(:,7))+ymax/10,'react.','HorizontalAlignment','center');
set(gca,'XTick',[0,2,5,10],'XColor','k','YColor','k');

title('DNA metylation');
xlabel('days after KRAB release');
ylabel('median DNAme density within amplicon');

saveas(fh,[parDIR,'medianDNAme_KRABtimecourse_release.pdf']);

%% avg DNAme at each position; rep combined
% actual index is like 50 -> 36 (-14 at the middle) -> 40 (+4 for spacing)
xaxlab = [10,20,50,60];
xax    = [10,20,40,50];
SizedFig(15,100);
visidx1 = [26,28,20,22]; % HDAC4
visidx2 = [27,29,21,23]; % HDAC4
% visidx1 = [1,24];
% visidx2 = [1,25];
% visidx1 = [14,15,16,7]; % KRAB rel2,5,10,irr
% visidx2 = [17,18,19,10]; % KRAB rel2,5,10,irr
% visidx1 = [5,NaN,6,14,15,16,7,12]; % KRAB nd,5d,rel2,5,10,irr,react
% visidx2 = [8,NaN,9,17,18,19,10,13]; % KRAB nd,5d,rel2,5,10,irr,react
visidx1 = [5,NaN,6,14,15,16,7,12]; % KRAB nd,5d,rel2,5,10,irr,react
visidx2 = [8,NaN,9,17,18,19,10,13]; % KRAB nd,5d,rel2,5,10,irr,react
% visidx1 = [6,14,15,16,7]; % 5d,rel2,5,10,irr
% visidx2 = [9,17,18,19,10]; % 5d,rel2,5,10,irr
DNAmeEFP1or2 = nan(2,length(visidx1));
DNAmeEFP1 = nan(2,length(visidx1));
disp('---');
for ii = 1:length(visidx1)
    if ~isnan(visidx1(ii))
        subplot(length(visidx1),1,ii);
        tmp1 = dnamet{visidx1(ii)};
        tmp1(tmp1 == 2) = NaN; tmp1(tmp1 == -1) = NaN; 
        DNAmeEFP1or2(1,ii) = sum((tmp1(:,6) == 0) | (tmp1(:,9) == 0))/size(tmp1,1);
        DNAmeEFP1(1,ii) = sum((tmp1(:,6) == 0))/size(tmp1,1);
        tmp1 = 1-nanmean(tmp1,1);
        tmp2 = dnamet{visidx2(ii)};
        tmp2(tmp2 == 2) = NaN; tmp2(tmp2 == -1) = NaN; 
        DNAmeEFP1or2(2,ii) = sum((tmp2(:,6) == 0) | (tmp2(:,9) == 0))/size(tmp2,1);
        DNAmeEFP1(2,ii) = sum((tmp2(:,6) == 0))/size(tmp2,1);
        tmp2 = 1-nanmean(tmp2,1);
%         tmp2 = 0;
        tmp = (tmp1+tmp2)/2;
        bar(tmp,1,'EdgeColor','none','FaceColor',[0.9,0.2,0.2]);
        hold on;
%         plot([6,9,11],tmp([6,9,11]),'bo','MarkerFaceColor','none','MarkerSize',6,'LineWidth',1);
    %     colormap(cmap);
        ylim([0,1]);
        set(gca,'XTick',xax,'XTickLabel',xaxlab);
        set(gca,'XColor','k','YColor','k');
    end
end
%% plot average time course of EFP1or2 met
ymax = 1.2;
SizedFig(20,20);
hold on;
plot([5,5],[0,ymax],'k--','LineWidth',1,'Color',[0.9,0.2,0.2]);
DNAmeEFP1or2 = DNAmeEFP1or2(:,[1,3:end]);
DNAmeEFP1or2avg = mean(DNAmeEFP1or2,1);
plot([0,5,7,10,15],DNAmeEFP1or2avg(1:end-2),'r-','LineWidth',2,'Color',[0.9,0.2,0.2]);
plot([0,5,7,10,15],DNAmeEFP1or2(:,1:end-2),'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','w','Color',[0.9,0.2,0.2]);

plot([19,22],DNAmeEFP1or2(:,[6,7])','ro','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','w','Color',[0.9,0.2,0.2]);


text(19,max(DNAmeEFP1or2(:,6))+ymax/10,'irr.','HorizontalAlignment','center');
text(22,max(DNAmeEFP1or2(:,7))+ymax/10,'react.','HorizontalAlignment','center');
set(gca,'XTick',[0,5,7,10,15],'XColor','k','YColor','k');

title('DNA metylation');
xlabel('days');
ylabel('Fraction DNAme on EFP1or2');
%% plot average time course of EFP1or2 met
ymax = 1.2;
SizedFig(20,20);
hold on;
plot([5,5],[0,ymax],'k--','LineWidth',1,'Color',[0.9,0.2,0.2]);
DNAmeEFP1 = DNAmeEFP1(:,3:end);
DNAmeEFP1or2avg = mean(DNAmeEFP1,1);
plot([5,7,10,15],DNAmeEFP1or2avg(1:end-2),'r-','LineWidth',2,'Color',[0.9,0.2,0.2]);
plot([5,7,10,15],DNAmeEFP1(:,1:end-2),'ro','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','w','Color',[0.9,0.2,0.2]);

plot([19,22],DNAmeEFP1(:,[6,7])','ro','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','w','Color',[0.9,0.2,0.2]);


text(19,max(DNAmeEFP1(:,6))+ymax/10,'irr.','HorizontalAlignment','center');
text(22,max(DNAmeEFP1(:,7))+ymax/10,'react.','HorizontalAlignment','center');
set(gca,'XTick',[0,5,7,10,15],'XColor','k','YColor','k');

title('DNA metylation');
xlabel('days');
ylabel('Fraction DNAme on EFP1or2');
%% avg DNAme at each position; rep combined; superimposed
% actual index is like 50 -> 36 (-14 at the middle) -> 40 (+4 for spacing)
% visidx1 = [14,15,16,7];
% visidx2 = [17,18,19,10];
visidx1 = [16];
visidx2 = [19];
dmetmat = [];
for ii = 1:length(visidx1)
    tmp1 = dnamet{visidx1(ii)};
    tmp1(tmp1 == 2) = NaN; tmp1(tmp1 == -1) = NaN; 
    tmp1 = 1-nanmean(tmp1,1);
    tmp2 = dnamet{visidx2(ii)};
    tmp2(tmp2 == 2) = NaN; tmp2(tmp2 == -1) = NaN; tmp2 = 1-nanmean(tmp2,1);
    dmetmat = [dmetmat;(tmp1+tmp2)/2];
end

xaxlab = [10,20,50,60];
xax    = [10,20,40,50];
cmap = [0.8,0.1,0.1 ; 0.9,0.5,0.5 ; 1.0,0.7,0.7; 1.0,0.9,0.9];
SizedFig(20,10);
hold on;
for ii = 1:size(dmetmat,1)
    tmp = dmetmat(size(dmetmat,1)-(ii-1),:);
    bar(tmp,1,'EdgeColor','none','FaceColor',cmap(ii,:));
    tmp([1:5,7,8,10,12:end]) = NaN;
    bar(tmp,1,'EdgeColor','none','FaceColor',cmap(ii+1,:));
%     colormap(cmap);
    ylim([0,1]);
    set(gca,'XTick',xax,'XTickLabel',xaxlab);
    set(gca,'XColor','k','YColor','k');
    
end




%% distribution of DNAme fraction
MedDNAmet = [];
xaxlab = [0,0.5,1];
xax    = [0,0.5,1];
SizedFig(10,100);
visidx1 = [5,6,14,15,16,7,12];
visidx2 = [8,9,17,18,19,10,13];
titles = {'no dox','dox 5days','dox5 release 2','dox5 release 5',...
    'dox5 release 10','irreversible','reactivated'};
for ii = 1:length(visidx1)
    if ii < 2
        subplot(length(visidx1)+1,1,ii);
    else
        subplot(length(visidx1)+1,1,ii+1);
    end
    tmp1 = dnamet{visidx1(ii)};
    tmp1(tmp1 == 2) = NaN; tmp1(tmp1 == -1) = NaN; tmp1 = 1-nanmean(tmp1,2);
    tmp2 = dnamet{visidx2(ii)};
    tmp2(tmp2 == 2) = NaN; tmp2(tmp2 == -1) = NaN; tmp2 = 1-nanmean(tmp2,2);
    histogram([tmp1;tmp2],20,'BinLimits',[0,1],'normalization','probability',...
        'FaceColor',[0.9,0.2,0.2]);
    alpha(1.0);
    set(gca,'XTick',xax,'XTickLabel',xaxlab);
    set(gca,'XColor','k','YColor','k');
    box off;
%     title(titles{ii});
    ylabel('density');
end
xlabel('fraction');

%% ----------------------------------------------------------------
%% -------------- two runs were combined and sorted ---------------
%% ----------------------------------------------------------------
%% heatmap, xax adjusted
dname = {'DNMT3B no dox','DNMT3B dox 12 days','DNMT3B/A6 no dox','DNMT3B/A6 dox 12 days',...
    'KRAB no dox rep1','KRAB dox 5 days rep1','irreversible rep1',...
    'KRAB no dox rep2','KRAB dox 5 days rep2','irreversible rep2',...
    'Abby control','reactivated rep1','reactivated rep2'};
% actual index is like 50 -> 36 (-14 at the middle) -> 40 (+4 for spacing)
xaxlab = [10,20,50,60];
xax    = [10,20,40,50];
cmap = [1.0,1.0,1.0 ; 0.8, 0.8, 0.8; 0.0, 0.0, 0.0];
cmap = [1.0,1.0,1.0 ; 0.9, 0.9, 0.9; 0.0, 0.0, 0.0];
cmap = [1.0,1.0,1.0 ; 0.0, 0.0, 0.0 ; 0.9, 0.9, 0.9; 0.5,0,0];
SizedFig(15,100);
visidx = 1:13; % ALL
visidx = [1,2,5:10,12,13]; % except A6
visidx = [1,2,5,6,7,12]; % rep1
% visidx = [1,2,8,9,10,13]; % rep2
for ii = 1:length(visidx)
    subplot(length(visidx),1,ii);
%     if ii <= 9
%     subplot(length(visidx),1,ii);
%     else
%     subplot(length(visidx),1,ii+1);
%     end
    tmp = dnamet{visidx(ii)};
%     imagesc(tmp);
    tmpidx = dnamet{visidx(ii)};
    tmpidx(tmpidx == 2) = NaN;
    tmpidx(tmpidx == -1) = NaN;
    tmpidx = 1-nanmean(tmpidx,2);
    [~,sortidx] = sort(nanmean(tmpidx,2));
    imagesc(tmp(sortidx,:));
    colormap(cmap);
%     xlim([0,100]);
    caxis([-1,2]);
    title(dname{visidx(ii)});
    set(gca,'XTick',xax,'XTickLabel',xaxlab);
    set(gca,'XColor','k','YColor','k');
end


%% median of DNAmet fraction
MedDNAmet = [];
SizedFig(12,100);
% visidx = 1:13; % ALL
% visidx = [1,2,5:10,12,13]; % except A6
visidx = [1,2,5,6,7,12]; % rep1
% visidx = [1,2,8,9,10,13]; % rep2
for ii = 1:length(visidx)
    subplot(length(visidx),1,ii);
    tmp = dnamet{visidx(ii)};
    tmp(tmp == 2) = NaN;
    tmp(tmp == -1) = NaN;
    tmp = 1-nanmean(tmp,2);
    histogram(tmp,20,'BinLimits',[0,1],'normalization','probability');
    set(gca,'XTick',xax,'XTickLabel',xaxlab);
    set(gca,'XColor','k','YColor','k');
    title(dname{visidx(ii)});
    num2str(round(nanmedian(tmp),3))
%     title(num2str(round(nanmedian(tmp),3)));
%     MedDNAmet = [MedDNAmet,nanmedian(tmp)];
%     title(num2str(round(sum(tmp > 0.5)/length(tmp),3)));
%     MedDNAmet = [MedDNAmet,sum(tmp > 0.5)/length(tmp)];
end
% MedDNAmet = cat(1,MedDNAmet([5,6,7,12,1,2]),[MedDNAmet([8,9,10,13]),NaN,NaN]);
% DNAmetscore = nanmean(MedDNAmet,1);

% figure;
% plot(tmp,'bo');

%%
%%
%%
%%

%% import EM-seq results, heatmap met

% L19	B9-DNMT3B -dox
% L20	B9-DNMT3B dox 12 days
% L21	A6-DNMT3B no dox
% L22	A6-DNMT3B dox 12 days
% L23	B9-KRAB no dox b1
% L24	B9-KRAB dox 5 days b1
% L25	B9-KRAB irrversible b1
% L26	B9-KRAB no dox b2
% L27	B9-KRAB dox 5 days b2
% L28	B9-KRAB irrversible b2
% -- 1:unmet, 0:met -- %

cpgidx = [31  ,38  ,51  ,72  ,84  ,103 ,105 ,129 ,142 ,153 ,170 ,189 ,192 ,198 ,208 ,213 ,222 ,243 ,256 ,258 ,275 ,291 ,324 ,339 ,346 ,371 ,381 ,398 ,403 ,430 ,439 ,442 ,444 ,448 ,466 ,468 ,477 ,486 ,527 ,530 ,563 ,590 ,603 ,605 ,609 ,612 ,615 ,622 ,626 ,634 ,644 ,647 ,652 ,661 ,665 ,667 ,674 ,681 ,685 ,706 ,729 ,731 ,734 ,737 ,744 ,749 ,758]+1;

for ii = 1:length(dnames)
    tmp = readmatrix([DIR2,dnames{ii}] ,'FileType', 'text');
    tmp = tmp(2:end,2:end);
    disp(size(tmp,1));
    tmp = tmp(:,cpgidx([1:26,27,28,29,30,(end-26):end]));
    tmp(:,[27,28,29,30]) = 2;
    if ii < 11
        jj = ii;
    else
        jj = ii + 1;
    end
    dnamet{jj} = cat(1,dnamet{jj},tmp);
end
%% ----------
%% bar plots
%% ----------
%% krab 5 days
idx1 = 5;
idx2 = 6;
xaxlab = [10,20,50,60];
xaxpos    = [10,20,40,50];
xax = 1:size(dnamet,2);
SizedFig(30,12);
hold on;
bar(xax,(dnamet(idx2,xax)+dnamet(idx2+3,xax))/2,1,'EdgeAlpha',0,'FaceAlpha',1.0);
bar(xax,(dnamet(idx1,xax)+dnamet(idx1+3,xax))/2,1,'EdgeAlpha',0,'FaceAlpha',1.0);
plot([xax;xax],[dnamet(idx2,xax);dnamet(idx2+3,xax)],'b-');
plot([xax;xax],[dnamet(idx1,xax);dnamet(idx1+3,xax)],'r-');
box on;
set(gca,'XTick',xaxpos,'XTickLabel',xaxlab);
set(gca,'XColor','k','YColor','k');
title('KRAB 5 days');
ylim([0,100]);
xlim([0.5,57.5]);
%% irrversible
idx1 = 5;
idx2 = 7;
xaxlab = [10,20,50,60];
xaxpos    = [10,20,40,50];
xax = 1:size(dnamet,2);
SizedFig(30,12);
hold on;
bar(xax,(dnamet(idx2,xax)+dnamet(idx2+0,xax))/2,1,'EdgeAlpha',0,'FaceAlpha',1.0);
bar(xax,(dnamet(idx1,xax)+dnamet(idx1+3,xax))/2,1,'EdgeAlpha',0,'FaceAlpha',1.0);
plot([xax;xax],[dnamet(idx2,xax);dnamet(idx2+0,xax)],'b-');
plot([xax;xax],[dnamet(idx1,xax);dnamet(idx1+3,xax)],'r-');
ylim([0,100]);
xlim([0.5,57.5]);
box on;
set(gca,'XTick',xaxpos,'XTickLabel',xaxlab);
set(gca,'XColor','k','YColor','k');
title('irreversible');
%% DNMT3B control
idx1 = 1;
idx2 = 2;
xaxlab = [10,20,50,60];
xaxpos    = [10,20,40,50];
xax = 1:size(dnamet,2);
SizedFig(30,12);
hold on;
bar(xax,dnamet(idx2,xax),1,'EdgeAlpha',0,'FaceAlpha',1.0);
bar(xax,dnamet(idx1,xax),1,'EdgeAlpha',0,'FaceAlpha',1.0);
ylim([0,100]);
xlim([0.5,57.5]);
box on;
set(gca,'XTick',xaxpos,'XTickLabel',xaxlab);
set(gca,'XColor','k','YColor','k');
    title('DNMT3B');
%% import EM-seq results, histogram of % met

% L19	B9-DNMT3B -dox
% L20	B9-DNMT3B dox 12 days
% L21	A6-DNMT3B no dox
% L22	A6-DNMT3B dox 12 days
% L23	B9-KRAB no dox b1
% L24	B9-KRAB dox 5 days b1
% L25	B9-KRAB irrversible b1
% L26	B9-KRAB no dox b2
% L27	B9-KRAB dox 5 days b2
% L28	B9-KRAB irrversible b2
% -- 1:unmet, 0:met -- %

cpgidx = [31  ,38  ,51  ,72  ,84  ,103 ,105 ,129 ,142 ,153 ,170 ,189 ,192 ,198 ,208 ,213 ,222 ,243 ,256 ,258 ,275 ,291 ,324 ,339 ,346 ,371 ,381 ,398 ,403 ,430 ,439 ,442 ,444 ,448 ,466 ,468 ,477 ,486 ,527 ,530 ,563 ,590 ,603 ,605 ,609 ,612 ,615 ,622 ,626 ,634 ,644 ,647 ,652 ,661 ,665 ,667 ,674 ,681 ,685 ,706 ,729 ,731 ,734 ,737 ,744 ,749 ,758]+1;

dnamet = {};
for ii = 1:length(dnames)
    tmp = readmatrix([DIR,dnames{ii}] ,'FileType', 'text');
    tmp = tmp(2:end,2:end);
    disp(size(tmp,1));
    tmp(tmp == -1) = NaN;
    tmp = tmp(:,cpgidx([1:26,end-26:end]));
    tmp = nanmean(tmp(:,1:26),2);
    tmp = (1 - tmp)*100; % % methylated
    dnamet = [dnamet,tmp];
end
%%
SizedFig(10,50);
% visidx = [1,2,8,9,10];
visidx = 1:13;
for ii = 1:length(visidx)
    subplot(5,1,ii);
    tmp = dnamet{visidx(ii)};
    disp((sum(tmp < 50) / length(tmp))*100);
    histogram(tmp,20,'BinLimits',[0,100],'Normalization','probability','EdgeColor','none','FaceColor','#0000ff');
    alpha(1.0);
    xlim([0,100]);
    set(gca,'XColor','k','YColor','k');
end


%% heatmap, xax adjusted, rep combined
dname = {'DNMT3B no dox','DNMT3B dox 12 days','DNMT3B/A6 no dox','DNMT3B/A6 dox 12 days',...
    'KRAB no dox rep1','KRAB dox 5 days rep1','irreversible rep1',...
    'KRAB no dox rep2','KRAB dox 5 days rep2','irreversible rep2',...
    'Abby control','reactivated rep1','reactivated rep2',...
    'KRAB dox5 release2 rep1','KRAB dox5 release5 rep1','KRAB dox5 release10 rep1',...
    'KRAB dox5 release2 rep2','KRAB dox5 release5 rep2','KRAB dox5 release10 rep2',...
    'HDAC4 reactivated rep1','HDAC4 reactivated rep2','HDAC4 irreversible rep1','HDAC4 irreversible rep2',...
    'DNMT3B irreversible rep1','DNMT3B irreversible rep2',...
    'HDAC4 no dox rep1','HDAC4 no dox rep2','HDAC4 dox 5days rep1','HDAC4 dox 5days rep2'};
% actual index is like 50 -> 36 (-14 at the middle) -> 40 (+4 for spacing)
xaxlab = [10,20,50,60];
xax    = [10,20,40,50];
cmap = [1.0,1.0,1.0 ; 0.0, 0.0, 0.0 ; 0.9, 0.9, 0.9; 0.5,0,0];
SizedFig(15,100);
visidx = [6,9; 14,17; 15,18; 16,19; 7,10]; % KRAB irreversible rep1&2
for jj = 1:size(visidx,1)
    subplot(size(visidx,1),1,jj);
    tmp = [];
    tmpidx = [];
    for ii = 1:size(visidx,2)
        tmp = [tmp;dnamet{visidx(jj,ii)}];
        tmpidx = [tmpidx;dnamet{visidx(jj,ii)}];
    end
    tmpidx(tmpidx == 2) = NaN;
    tmpidx(tmpidx == -1) = NaN;
    tmpidx = 1-nanmean(tmpidx,2);
    [~,sortidx] = sort(nanmean(tmpidx,2));
    imagesc(tmp(sortidx,:));
    colormap(cmap);
    %     xlim([0,100]);
    caxis([-1,2]);
    title(dname{visidx(jj,ii)});
    set(gca,'XTick',xax,'XTickLabel',xaxlab);
    set(gca,'XColor','k','YColor','k');
end

%%
%%
%% call C of CpG and G of CpG separately
%%
%% import EM-seq results, heatmap met

% -- 1:unmet, 0:met -- %

cpgidx = [31  ,38  ,51  ,72  ,84  ,103 ,105 ,129 ,142 ,153 ,170 ,189 ,192 ,198 ,208 ,213 ,222 ,243 ,256 ,258 ,275 ,291 ,324 ,339 ,346 ,371 ,381 ,398 ,403 ,430 ,439 ,442 ,444 ,448 ,466 ,468 ,477 ,486 ,527 ,530 ,563 ,590 ,603 ,605 ,609 ,612 ,615 ,622 ,626 ,634 ,644 ,647 ,652 ,661 ,665 ,667 ,674 ,681 ,685 ,706 ,729 ,731 ,734 ,737 ,744 ,749 ,758];

dnametC = {};
for ii = 1:length(dir1idx)
    tmp1 = [];
    if ~isnan(dir1idx(ii))
        cpgidx_tmp = cpgidx + 1;
        tmp1 = readmatrix([DIR1,dnames{dir1idx(ii)}] ,'FileType', 'text');
        tmp1 = tmp1(2:end,2:end);
        tmp1 = tmp1(:,cpgidx_tmp([1:26,27,28,29,30,(end-26):end]));
        tmp1(:,[27,28,29,30]) = 2;
    end
    tmp2 = [];
    if ~isnan(dir2idx(ii))
        cpgidx_tmp = cpgidx;
        tmp2 = readmatrix([DIR2,dnames2{dir2idx(ii)}] ,'FileType', 'text');
        tmp2 = tmp2(2:end,2:end);
        tmp2 = tmp2(:,cpgidx_tmp([1:26,27,28,29,30,(end-26):end]));
        tmp2(:,[27,28,29,30]) = 2;
    end
    tmp3 = [];
    if ~isnan(dir3idx(ii))
        cpgidx_tmp = cpgidx;
        tmp2 = readmatrix([DIR3,dnames3{dir3idx(ii)}] ,'FileType', 'text');
        tmp2 = tmp2(2:end,2:end);
        tmp2 = tmp2(:,cpgidx_tmp([1:26,27,28,29,30,(end-26):end]));
        tmp2(:,[27,28,29,30]) = 2;
    end
    tmp = cat(1,tmp1,tmp2,tmp3);
    disp([size(tmp1,1),size(tmp2,1),size(tmp3,1),size(tmp,1)]);
    dnametC = [dnametC,tmp];
end


dnametG = {};
for ii = 1:length(dir1idx)
    tmp1 = [];
    if ~isnan(dir1idx(ii))
        cpgidx_tmp = cpgidx + 2;
        tmp1 = readmatrix([DIR1,dnames{dir1idx(ii)}] ,'FileType', 'text');
        tmp1 = tmp1(2:end,2:end);
        tmp1 = tmp1(:,cpgidx_tmp([1:26,27,28,29,30,(end-26):end]));
        tmp1(:,[27,28,29,30]) = 2;
    end
    tmp2 = [];
    if ~isnan(dir2idx(ii))
        cpgidx_tmp = cpgidx + 1;
        tmp2 = readmatrix([DIR2,dnames2{dir2idx(ii)}] ,'FileType', 'text');
        tmp2 = tmp2(2:end,2:end);
        tmp2 = tmp2(:,cpgidx_tmp([1:26,27,28,29,30,(end-26):end]));
        tmp2(:,[27,28,29,30]) = 2;
    end
    tmp3 = [];
    if ~isnan(dir3idx(ii))
        cpgidx_tmp = cpgidx + 1;
        tmp3 = readmatrix([DIR3,dnames3{dir3idx(ii)}] ,'FileType', 'text');
        tmp3 = tmp3(2:end,2:end);
        tmp3 = tmp3(:,cpgidx_tmp([1:26,27,28,29,30,(end-26):end]));
        tmp3(:,[27,28,29,30]) = 2;
    end
    tmp = cat(1,tmp1,tmp2,tmp3);
    disp([size(tmp1,1),size(tmp2,1),size(tmp3,1),size(tmp,1)]);
    dnametG = [dnametG,tmp];
end

%% heatmap, xax adjusted
% actual index is like 50 -> 36 (-14 at the middle) -> 40 (+4 for spacing)
xaxlab = [10,20,50,60];
xax    = [10,20,40,50];
cmap = [1.0,1.0,1.0 ; 0.0, 0.0, 0.0 ; 0.9, 0.9, 0.9; 0.5,0,0];
SizedFig(15,100);
visidx = 1:29;
% visidx = [5,6,14,15,16,7,12]; % KRAB recruitment rep1
% visidx = [8,9,17,18,19,10,13]; % KRAB recruitment rep2
for ii = 1:length(visidx)
    subplot(length(visidx),1,ii);
%     if ii <= 9
%     subplot(length(visidx),1,ii);
%     else
%     subplot(length(visidx),1,ii+1);
%     end
    tmp = dnametC{visidx(ii)};
    disp(sum(sum(dnametC{visidx(ii)}-dnametG{visidx(ii)} ~= 0)));
    imagesc(tmp);
    colormap(cmap);
%     xlim([0,100]);
    caxis([-1,2]);
    title(dname{visidx(ii)});
    set(gca,'XTick',xax,'XTickLabel',xaxlab);
    set(gca,'XColor','k','YColor','k');
end
