%% params
DIR = 'D:\Taihei\CUTnRUN\20220927_KRAB_CNR_H3K9me3\dedup_scaled_chr19\';
fnames = dir(DIR);
fnames = fnames(3:end);

%% indices 5kb
downstreamend = 55115765; % upstreamstart = 55115767;
reporterstart = 55115766;
reporterend = 55120420;

downstream = (downstreamend-5000*14):5000:downstreamend;
downstream = [downstream'-4999,downstream'];

upstream = (reporterend+1):5000:(reporterend+5000*15);
upstream = [upstream',upstream'+4999];

stedidx = [downstream;[reporterstart,reporterend];upstream];
%% indices 1kb
downstreamend = 55115765; % upstreamstart = 55115767;
reporterstart = 55115766;
reporterend = 55120420;
stepsize = 1000;
tmp = (reporterstart:931:reporterend)';
tmp = [tmp,[tmp(2:end);reporterend]];

downstream = (downstreamend-stepsize*(15*5-1)):stepsize:downstreamend;
downstream = [downstream'-(stepsize-1),downstream'];

upstream = (reporterend+1):stepsize:(reporterend+stepsize*15*5);
upstream = [upstream',upstream'+(stepsize-1)];

stedidx = [downstream;tmp;upstream];
%% indices 100bb
downstreamend = 55115765; % upstreamstart = 55115767;
reporterstart = 55115766;
reporterend = 55120420;
stepsize = 100;
tmp = (reporterstart:99:(reporterend-5))';
tmp = [tmp,[tmp(2:end);reporterend]];

downstream = (downstreamend-stepsize*(100-1)):stepsize:downstreamend;
downstream = [downstream'-(stepsize-1),downstream'];

upstream = (reporterend+1):stepsize:(reporterend+stepsize*100);
upstream = [upstream',upstream'+(stepsize-1)];

stedidx = [downstream;tmp;upstream];
%% import files, calculate local normalized reads
MarkProf = zeros(size(stedidx,1),length(fnames));
for ii = 1:length(fnames)
    tst = readmatrix([DIR,fnames(ii).name]);
    for jj = 1:size(stedidx,1)
        idx = tst(:,2) >=stedidx(jj,1) & tst(:,3) <= stedidx(jj,2);
%         MarkProf(jj,ii) = sum(tst(idx,4));
        MarkProf(jj,ii) = sum(tst(idx,4))/((stedidx(jj,2) - stedidx(jj,1)+1)/1000); % /kb normalization
    end
    
end

permuteidx = [12,5,7,15,9,3,1, 13,6,8,16,10,4,2];
MarkProf = MarkProf(:,permuteidx);
for ii = 1:size(permuteidx,2)
    disp([num2str(ii),':',fnames(permuteidx(ii)).name]);
end
MarkProfAvg = (MarkProf(:,1:7)+MarkProf(:,8:14))/2;
% MarkProfH3K9meSd = std(cat(3,MarkProfH3K9me(:,1:7),MarkProfH3K9me(:,8:14)));
%%
%% ----------------------------------------------------------------------------------------
%% --------------------------------         plots          --------------------------------
%% ----------------------------------------------------------------------------------------
%% line plots
xax = -75:5:75;
 xax = -77:1:77;
 figure;
subplot(3,1,1);
plot(repmat(xax',[1 7]),MarkProfH3K9me(:,1:7),'o-');
subplot(3,1,2);
plot(repmat(xax',[1 7]),MarkProfH3K9me(:,8:14),'o-');
subplot(3,1,3);
plot(xax,MarkProf(:,11),'o-');

%% 
xax = -75:5:75;
%  xax = -77:1:77;
%  xax = (-123:1:123)*0.1; % xax for 1kb

figure;
% subplot(3,1,1);
for ii = 1:7
    subplot(7,1,ii);
%     bar(xax',flipud(MarkProfH3K9meAvg(:,ii)),1,'EdgeAlpha',0,'FaceAlpha',1);
    bar(xax',flipud(MarkProfAvg(:,ii)),1,'EdgeAlpha',0,'FaceAlpha',1);
    ylim([0,10]);
    xlim([-32.5,32.5]);
%     xlim([-14,14]);
end
% subplot(3,1,2);
% plot(repmat(xax',[1 7]),MarkProfH3K9me(:,8:14),'o-');
% subplot(3,1,3);
% plot(xax,MarkProf(:,11),'o-');