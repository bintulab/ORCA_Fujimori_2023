%%
pardir = 'D:\Taihei\CUTnRUN\20230720_KRAB_CNR_H3K9me2_H3K27me3\';
permuteidx =  [12,3,7,17,9,1,5,15];
permuteidx =  [permuteidx,permuteidx+1];
%% GAPDH 
% params
DIR = [pardir,'dedup_scaled_chr12\'];
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
MarkProfGAPDH = MarkProf(:,permuteidx);

%% RPL10A 
% params
DIR = [pardir,'dedup_scaled_chr6\'];
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
MarkProfRPL10A = MarkProf(:,permuteidx);

%% NXF2B chrX:102,360,395-102,440,008
% params
DIR = [pardir,'dedup_scaled_chrX\'];
fnames = dir(DIR);
fnames = fnames(3:end);

% indices
stedidx = [102360395,102440008];

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
MarkProfNXF2B = MarkProf(:,permuteidx);

%% chr21 DSCR4 hg38 21:38,054,011-38,121,360
% params
DIR = [pardir,'dedup_scaled_chr21\'];
fnames = dir(DIR);
fnames = fnames(3:end);

% indices

stedidx = [38054011,38121360];

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
MarkProfDSCR4 = MarkProf(:,permuteidx);

%% reporter 
% params
DIR = [pardir,'dedup_scaled_chr19\'];
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
MarkProfRep = MarkProf(:,permuteidx);
%%
%% ----------------------------------------------------------------------------------------
%% --------------------------------         plots          --------------------------------
%% ----------------------------------------------------------------------------------------
MarkProfGAPDHAvg = (MarkProfGAPDH(:,1:7)+MarkProfGAPDH(:,8:14))/2;
MarkProfGAPDHSd = std(cat(3,MarkProfGAPDH(:,1:7),MarkProfGAPDH(:,8:14)));
%% bar plots
SizedFig(50,20);
visidx = [1,2,3,4]; barpos = 2;
% visidx = [5,6,7,8]; barpos = 5;
range = 1:length(visidx);
xax = -1;
hold on;
repoffset = 8;

for ii = 1:5
    switch ii 
        case 1
            tmp = MarkProfGAPDH;
            gname = 'GAPDH';
        case 2
            tmp = MarkProfRPL10A;
            gname = 'RPL10A';
        case 3
            tmp = MarkProfDSCR4;
            gname = 'DSCR4';
        case 4
            tmp = MarkProfNXF2B;
            gname = 'NXF2B';
        case 5
            tmp = MarkProfRep;
            gname = 'reporter +/- 30kb';
    end
    offset = xax(end)+1; 
    xax = range+offset;
    bar(xax, (tmp(visidx)+tmp(visidx+repoffset))/2);
    plot(xax, tmp(visidx),'ko','MarkerFaceColor','w');
    plot(xax, tmp(visidx+repoffset),'ko','MarkerFaceColor','w');
    plot([xax(1)-0.5,xax(end)+0.5],[barpos,barpos],'k-');
    text(mean(xax),barpos*1.1,gname,'HorizontalAlignment','center','FontSize',15);
end
set(gca,'XTick',[],'XColor','k','YColor','k','FontSize',15);
% set(gca,'Yscale','log');
ylim([0,barpos*1.2]);
%%
tmp = MarkProfMEI4;
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
    bar(xax',MarkProfGAPDHAvg(:,ii)-MarkProfGAPDHAvg(:,1),1,'EdgeAlpha',0,'FaceAlpha',1);
    ylim([-1,40]);
    xlim([-32.5,32.5]);
end
% subplot(3,1,2);
% plot(repmat(xax',[1 7]),MarkProf(:,8:14),'o-');
% subplot(3,1,3);
% plot(xax,MarkProf(:,11),'o-');