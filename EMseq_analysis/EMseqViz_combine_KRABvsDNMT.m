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
fnames = fnames(3:end);
fnames = fnames(11:12);
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
dname = {'KRAB no dox rep1','KRAB no dox rep2',...
    'KRAB irreversible rep1','KRAB irreversible rep2',...
    'KRAB reactivated rep1','KRAB reactivated rep2',...
    'DNMT3B no dox rep1','DNMT3B no dox rep2',...
    'DNMT3B irreversible rep1','DNMT3B irreversible rep1'};
dir1idx = [  5,  8,  7, 10,NaN,NaN,  1,NaN,NaN,NaN];
dir2idx = [  5,  8,  7,NaN, 11, 12,NaN,NaN,NaN,NaN];
dir3idx = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,  1,  2];
%% import EM-seq results, heatmap met

% -- 1:unmet, 0:met -- %

cpgidx = [31  ,38  ,51  ,72  ,84  ,103 ,105 ,129 ,142 ,153 ,170 ,189 ,192 ,198 ,208 ,213 ,222 ,243 ,256 ,258 ,275 ,291 ,324 ,339 ,346 ,371 ,381 ,398 ,403 ,430 ,439 ,442 ,444 ,448 ,466 ,468 ,477 ,486 ,527 ,530 ,563 ,590 ,603 ,605 ,609 ,612 ,615 ,622 ,626 ,634 ,644 ,647 ,652 ,661 ,665 ,667 ,674 ,681 ,685 ,706 ,729 ,731 ,734 ,737 ,744 ,749 ,758]+1;

dnamet = cell(length(dir1idx),1);
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
        tmp2 = readmatrix([DIR3,dnames3{dir3idx(ii)}] ,'FileType', 'text');
        tmp2 = tmp2(2:end,2:end);
        tmp2 = tmp2(:,cpgidx([1:26,27,28,29,30,(end-26):end]));
        tmp2(:,[27,28,29,30]) = 2;
    end
    tmp = cat(1,tmp1,tmp2,tmp3);
    disp([size(tmp1,1),size(tmp2,1),size(tmp3,1),size(tmp,1)]);
    dnamet{ii} = tmp;
end

%% calc median of DNAmet fraction
MedDNAmet = [];
SizedFig(15,100);
visidx = 1:10;
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
    MedDNAmet = [MedDNAmet,nanmedian(tmp)];
end
MedDNAmet = cat(1,MedDNAmet(1:2:end),MedDNAmet(2:2:end));
DNAmetscore = nanmean(MedDNAmet,1);

% figure;
% plot(tmp,'bo');

%%  average median DNAme
ymax = 0.75;
SizedFig(20,30);
hold on;
cond_names = {'no dox','irreversible','reactivated',...
    'no dox','irreversible'};
xaxtmp = [0,1,2,3.5,4.5];
bar(xaxtmp,DNAmetscore(1:end),'FaceColor',[0.8,0.8,0.8]);
plot(xaxtmp,MedDNAmet(:,1:end),'o','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','w','Color',[0.9,0.2,0.2]);

set(gca,'XTick',[0,5,7,10,15],'XColor','k','YColor','k');

ylabel('median DNAme fraction');
set(gca,'XTick',xaxtmp(1:end),'XTickLabel',cond_names);
xtickangle(-45);
ylim([0,1.2]);
%% avg DNAme at each position; rep combined
% actual index is like 50 -> 36 (-14 at the middle) -> 40 (+4 for spacing)
xaxlab = [10,20,50,60];
xax    = [10,20,40,50];
cmap = [1.0,1.0,1.0 ; 0.0, 0.0, 0.0 ; 0.9, 0.9, 0.9; 0.5,0,0];
SizedFig(15,100);
titles = {'KRAB no dox','KRAB irreversible','KRAB reactivated',...
    'DNMT3B no dox','DNMT3B irreversible'};
visidx1 = 1:2:10; % rep1
visidx2 = 2:2:10; % rep2
for ii = 1:length(visidx1)
    subplot(length(visidx1),1,ii);
    tmp1 = dnamet{visidx1(ii)};
    tmp1(tmp1 == 2) = NaN; tmp1(tmp1 == -1) = NaN; tmp1 = 1-nanmean(tmp1,1);
    tmp2 = dnamet{visidx2(ii)};
    tmp2(tmp2 == 2) = NaN; tmp2(tmp2 == -1) = NaN; tmp2 = 1-nanmean(tmp2,1);
    if isempty(tmp2)
        tmp2 = tmp1;
    end
    bar((tmp1+tmp2)/2,1,'EdgeColor','none','FaceColor',[0.9,0.2,0.2]);
%     colormap(cmap);
    ylim([0,1]);
    title(titles{ii});
    set(gca,'XTick',xax,'XTickLabel',xaxlab);
    set(gca,'XColor','k','YColor','k');
end
%% median of DNAmet fraction
MedDNAmet = [];
xaxlab = [0,0.5,1];
xax    = [0,0.5,1];
SizedFig(15,100);
titles = {'KRAB no dox','KRAB irreversible','KRAB reactivated',...
    'DNMT3B no dox','DNMT3B irreversible'};
visidx1 = 1:2:10; % rep1
visidx2 = 2:2:10; % rep2
for ii = 1:length(visidx1)
    subplot(length(visidx1),1,ii);
    tmp1 = dnamet{visidx1(ii)};
    tmp1(tmp1 == 2) = NaN; tmp1(tmp1 == -1) = NaN; tmp1 = 1-nanmean(tmp1,2);
    tmp2 = dnamet{visidx2(ii)};
    tmp2(tmp2 == 2) = NaN; tmp2(tmp2 == -1) = NaN; tmp2 = 1-nanmean(tmp2,2);
    histogram([tmp1;tmp2],20,'BinLimits',[0,1],'normalization','probability',...
        'FaceColor',[0.9,0.2,0.2]);
    alpha(1.0);
    set(gca,'XTick',xax,'XTickLabel',xaxlab);
    set(gca,'XColor','k','YColor','k');
    title(titles{ii});
end
