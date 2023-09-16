%% GAPDH 
% params
DIR = 'D:\Taihei\CUTnRUN\20220927_KRAB_CNR_H3K9me3\dedup_scaled_chr12\';
fnames = dir(DIR);
fnames = fnames(3:end);

% indices
stedidx = [6534754,6538339];

% import files, calculate local normalized reads
MarkProf = zeros(size(stedidx,1),length(fnames));
for ii = 1:length(fnames)
    tst = readmatrix([DIR,fnames(ii).name]);
    for jj = 1:size(stedidx,1)
        idx = tst(:,2) >=stedidx(jj,1) & tst(:,3) <= stedidx(jj,2);
        %         MarkProf(jj,ii) = sum(tst(idx,4));
        MarkProf(jj,ii) = sum(tst(idx,4))/((stedidx(jj,2) - stedidx(jj,1)+1)/1000); % /kb normalization
    end
    disp(ii);
end

permuteidx = [12,5,7,15,9,3,1, 13,6,8,16,10,4,2];
MarkProfGAPDH = MarkProf(:,permuteidx);

%% RPL10A 
% params
DIR = 'D:\Taihei\CUTnRUN\20220927_KRAB_CNR_H3K9me3\dedup_scaled_chr6\';
fnames = dir(DIR);
fnames = fnames(3:end);

% indices
stedidx = [35468401,35470780];

% import files, calculate local normalized reads
MarkProf = zeros(size(stedidx,1),length(fnames));
for ii = 1:length(fnames)
    tst = readmatrix([DIR,fnames(ii).name]);
    for jj = 1:size(stedidx,1)
        idx = tst(:,2) >=stedidx(jj,1) & tst(:,3) <= stedidx(jj,2);
        %         MarkProf(jj,ii) = sum(tst(idx,4));
        MarkProf(jj,ii) = sum(tst(idx,4))/((stedidx(jj,2) - stedidx(jj,1)+1)/1000); % /kb normalization
    end
    disp(ii);
end

permuteidx = [12,5,7,15,9,3,1, 13,6,8,16,10,4,2];
MarkProfRPL10A = MarkProf(:,permuteidx);

%% CXCL12 
% params
DIR = 'D:\Taihei\CUTnRUN\20220927_KRAB_CNR_H3K9me3\dedup_scaled_chr10\';
fnames = dir(DIR);
fnames = fnames(3:end);

% indices
stedidx = [44372888,44385097];

% import files, calculate local normalized reads
MarkProf = zeros(size(stedidx,1),length(fnames));
for ii = 1:length(fnames)
    tst = readmatrix([DIR,fnames(ii).name]);
    for jj = 1:size(stedidx,1)
        idx = tst(:,2) >=stedidx(jj,1) & tst(:,3) <= stedidx(jj,2);
        %         MarkProf(jj,ii) = sum(tst(idx,4));
        MarkProf(jj,ii) = sum(tst(idx,4))/((stedidx(jj,2) - stedidx(jj,1)+1)/1000); % /kb normalization
    end
    disp(ii);
end

permuteidx = [12,5,7,15,9,3,1, 13,6,8,16,10,4,2];
MarkProfCXCL12 = MarkProf(:,permuteidx);


%% ZNF14 
% params
DIR = 'D:\Taihei\CUTnRUN\20220927_KRAB_CNR_H3K9me3\dedup_scaled_chr19\';
fnames = dir(DIR);
fnames = fnames(3:end);

% indices
stedidx = [19710472,19733112];

% import files, calculate local normalized reads
MarkProf = zeros(size(stedidx,1),length(fnames));
for ii = 1:length(fnames)
    tst = readmatrix([DIR,fnames(ii).name]);
    for jj = 1:size(stedidx,1)
        idx = tst(:,2) >=stedidx(jj,1) & tst(:,3) <= stedidx(jj,2);
        %         MarkProf(jj,ii) = sum(tst(idx,4));
        MarkProf(jj,ii) = sum(tst(idx,4))/((stedidx(jj,2) - stedidx(jj,1)+1)/1000); % /kb normalization
    end
    disp(ii);
end

permuteidx = [12,5,7,15,9,3,1, 13,6,8,16,10,4,2];
MarkProfZNF14 = MarkProf(:,permuteidx);

%% reporter 
% params
DIR = 'D:\Taihei\CUTnRUN\20220927_KRAB_CNR_H3K9me3\dedup_scaled_chr19\';
fnames = dir(DIR);
fnames = fnames(3:end);
% indices
downstreamend = 55115765; % upstreamstart = 55115767;
reporterstart = 55115766;
reporterend = 55120420;
segnum = 6;
downstream = (downstreamend-5000*(segnum-1)):5000:downstreamend;
downstream = [downstream'-4999,downstream'];
upstream = (reporterend+1):5000:(reporterend+5000*segnum);
upstream = [upstream',upstream'+4999];
stedidx = [downstream(1),upstream(end)];

% import files, calculate local normalized reads
MarkProf = zeros(size(stedidx,1),length(fnames));
for ii = 1:length(fnames)
    tst = readmatrix([DIR,fnames(ii).name]);
    for jj = 1:size(stedidx,1)
        idx = tst(:,2) >=stedidx(jj,1) & tst(:,3) <= stedidx(jj,2);
        %         MarkProf(jj,ii) = sum(tst(idx,4));
        MarkProf(jj,ii) = sum(tst(idx,4))/((stedidx(jj,2) - stedidx(jj,1)+1)/1000); % /kb normalization
    end
    disp(ii);
end

permuteidx = [12,5,7,15,9,3,1, 13,6,8,16,10,4,2];
MarkProfRep = MarkProf(:,permuteidx);
%%
%% ----------------------------------------------------------------------------------------
%% --------------------------------         plots          --------------------------------
%% ----------------------------------------------------------------------------------------
MarkProfH3K9meGAPDHAvg = (MarkProfH3K9meGAPDH(:,1:7)+MarkProfH3K9meGAPDH(:,8:14))/2;
MarkProfH3K9meGAPDHSd = std(cat(3,MarkProfH3K9meGAPDH(:,1:7),MarkProfH3K9meGAPDH(:,8:14)));
%% bar plots
SizedFig(50,20);
visidx = [1,3,4,5];
% visidx = 1:7;
range = 1:length(visidx);
xax = -1;
hold on;

for ii = 1:5
    switch ii 
        case 1
            tmp = MarkProfGAPDH;
            gname = 'GAPDH';
        case 2
            tmp = MarkProfRPL10A;
            gname = 'RPL10A';
        case 3
            tmp = MarkProfCXCL12;
            gname = 'CXCL12';
        case 4
            tmp = MarkProfZNF14;
            gname = 'ZNF14';
        case 5
            tmp = MarkProfRep;
            gname = 'reporter +/- 30kb';
    end
    offset = xax(end)+1; 
    xax = range+offset;
    bar(xax, (tmp(visidx)+tmp(visidx+7))/2);
    plot(xax, tmp(visidx),'ko','MarkerFaceColor','w');
    plot(xax, tmp(visidx+7),'ko','MarkerFaceColor','w');
    plot([xax(1)-0.5,xax(end)+0.5],[300,300],'k-');
    text(mean(xax),340,gname,'HorizontalAlignment','center','FontSize',15);
end
set(gca,'XTick',[],'XColor','k','YColor','k','FontSize',15);
ylim([0,10]);
%%
tmp = MarkProfH3K9meMEI4;
offset = 26; 
xax = range+offset;
bar(xax, (tmp([1,3,4,5])+tmp([1,3,4,5]+7))/2);
plot(xax, tmp([1,3,4,5]),'bo');
plot(xax, tmp([1,3,4,5]+7),'bo');

%% 
xax = -75:5:75;
figure;
% subplot(3,1,1);
for ii = 1:7
    subplot(1,7,ii);
    bar(xax',MarkProfH3K9meGAPDHAvg(:,ii)-MarkProfH3K9meGAPDHAvg(:,1),1,'EdgeAlpha',0,'FaceAlpha',1);
    ylim([-1,40]);
    xlim([-32.5,32.5]);
end
% subplot(3,1,2);
% plot(repmat(xax',[1 7]),MarkProfH3K9me(:,8:14),'o-');
% subplot(3,1,3);
% plot(xax,MarkProf(:,11),'o-');

%% --------------------------------
%% FEZF2 
% params
DIR = 'D:\Taihei\CUTnRUN\dedup_scaled_chr3\';
fnames = dir(DIR);
fnames = fnames(3:end);

% indices
stedidx = [62369681,62374324];

% import files, calculate local normalized reads
MarkProf = zeros(size(stedidx,1),length(fnames));
for ii = 1:length(fnames)
    tst = readmatrix([DIR,fnames(ii).name]);
    for jj = 1:size(stedidx,1)
        idx = tst(:,2) >=stedidx(jj,1) & tst(:,3) <= stedidx(jj,2);
        MarkProf(jj,ii) = sum(tst(idx,4));
    end
    disp(ii);
end

permuteidx = [12,5,7,15,9,3,1, 13,6,8,16,10,4,2];
MarkProfH3K9meFEZF2 = MarkProf(:,permuteidx);

%% MEI4 
% params
DIR = 'D:\Taihei\CUTnRUN\dedup_scaled_chr6\';
fnames = dir(DIR);
fnames = fnames(3:end);

% indices
stedidx = [77650274,77927028];

% import files, calculate local normalized reads
MarkProf = zeros(size(stedidx,1),length(fnames));
for ii = 1:length(fnames)
    tst = readmatrix([DIR,fnames(ii).name]);
    for jj = 1:size(stedidx,1)
        idx = tst(:,2) >=stedidx(jj,1) & tst(:,3) <= stedidx(jj,2);
        MarkProf(jj,ii) = sum(tst(idx,4));
    end
    disp(ii);
end

permuteidx = [12,5,7,15,9,3,1, 13,6,8,16,10,4,2];
MarkProfH3K9meMEI4 = MarkProf(:,permuteidx);
