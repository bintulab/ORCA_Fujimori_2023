%% params
DIR = '/Users/tfuji/Documents/tempanalysis/EM-seq/20230728_EM-seq_B9-KRAB_5d_irr_byAbby/heatmaps/';
% DIR = '/Users/tfuji/Documents/tempanalysis/EM-seq/20230723_EM-seq_B9-KRAB_5d_irr/';
fnames = dir(DIR);
fnames = fnames(3:end);
dnames = {};
for ii = 1:length(fnames)
    tmpcell = split(fnames(ii).name,'.');
    if strcmp(tmpcell{end},'matrix') && strcmp(tmpcell{end-1},'clustered')
        dnames = [dnames,fnames(ii).name];
    end
end
%% import EM-seq results

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

% for ii = 3%:length(dnames)
%     tmp = readmatrix([DIR,dnames{ii}] ,'FileType', 'text');
%     tmp = tmp(2:end,2:end);
%     tmp(tmp == -1) = NaN;
%     tmp = nanmean(tmp,1);
%     tmp = (1 - tmp)*100; % % methylated
% end
% cpgidx = ~isnan(tmp);
cpgidx = [31  ,38  ,51  ,72  ,84  ,103 ,105 ,129 ,142 ,153 ,170 ,189 ,192 ,198 ,208 ,213 ,222 ,243 ,256 ,258 ,275 ,291 ,324 ,339 ,346 ,371 ,381 ,398 ,403 ,430 ,439 ,442 ,444 ,448 ,466 ,468 ,477 ,486 ,527 ,530 ,563 ,590 ,603 ,605 ,609 ,612 ,615 ,622 ,626 ,634 ,644 ,647 ,652 ,661 ,665 ,667 ,674 ,681 ,685 ,706 ,729 ,731 ,734 ,737 ,744 ,749 ,758]+1;

dnamet = [];
for ii = 1:length(dnames)
    tmp = readmatrix([DIR,dnames{ii}] ,'FileType', 'text');
    tmp = tmp(2:end,2:end);
    disp(size(tmp,1));
    tmp(tmp == -1) = NaN;
    tmp = nanmean(tmp,1);
    tmp = (1 - tmp)*100; % % methylated
    
%     tmp = tmp(cpgidx);
    tmp = tmp(:,cpgidx([1:26,27,28,29,30,(end-26):end]));
    tmp(:,[27,28,29,30]) = 0;

    dnamet = [dnamet;tmp];
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
bar(xax,(dnamet(idx2,xax)+dnamet(idx2+3,xax))/2,1,'EdgeAlpha',0,'FaceAlpha',1.0);
bar(xax,(dnamet(idx1,xax)+dnamet(idx1+3,xax))/2,1,'EdgeAlpha',0,'FaceAlpha',1.0);
plot([xax;xax],[dnamet(idx2,xax);dnamet(idx2+3,xax)],'b-');
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
visidx = [1,2,8,9,10];
% visidx = [1,2,5,6,7];
% visidx = [3,4,5,7];
for ii = 1:length(visidx)
    subplot(5,1,ii);
    tmp = dnamet{visidx(ii)};
    disp((sum(tmp < 50) / length(tmp))*100);
    histogram(tmp,20,'BinLimits',[0,100],'Normalization','probability','EdgeColor','none','FaceColor','#0000ff');
    alpha(1.0);
    xlim([0,100]);
    set(gca,'XColor','k','YColor','k');
end


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

dnamet = {};
for ii = 1:length(dnames)
    tmp = readmatrix([DIR,dnames{ii}] ,'FileType', 'text');
    tmp = tmp(2:end,2:end);
    disp(size(tmp,1));
%     tmp = tmp(:,[1:cpgidx(27),cpgidx(end-26):end]);
    tmp = tmp(:,cpgidx([1:26,27,28,29,30,(end-26):end]));
    tmp(:,[27,28,29,30]) = 2;
    dnamet = [dnamet,tmp];
end
%% heatmap, xax adjusted
% actual index is like 50 -> 36 (-14 at the middle) -> 40 (+4 for spacing)
xaxlab = [10,20,50,60];
xax    = [10,20,40,50];
cmap = [1.0,1.0,1.0 ; 0.8, 0.8, 0.8; 0.0, 0.0, 0.0];
cmap = [1.0,1.0,1.0 ; 0.9, 0.9, 0.9; 0.0, 0.0, 0.0];
cmap = [1.0,1.0,1.0 ; 0.0, 0.0, 0.0 ; 0.9, 0.9, 0.9; 0.5,0,0];
SizedFig(15,50);
visidx = [1,2,5,6,7];
% visidx = [1,2,8,9,10];
% visidx = [3,4,5,7];
for ii = 1:length(visidx)
    subplot(5,1,ii);
    tmp = dnamet{visidx(ii)};
    imagesc(tmp);
    colormap(cmap);
%     xlim([0,100]);
    caxis([-1,2]);
    set(gca,'XTick',xax,'XTickLabel',xaxlab);
    set(gca,'XColor','k','YColor','k');
end

