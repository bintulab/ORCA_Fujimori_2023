%% params
DIR = 'D:\Taihei\CUTnRUN\20230720_KRAB_CNR_H3K9me2_H3K27me3\dedup_scaled_chr19\';
fnames = dir(DIR);
fnames = fnames(3:end);
%%
for ii = 1:length(fnames)
    disp([num2str(ii),':',fnames(ii).name]);
end
permuteidx =  [12,3,7,17,9,1,5,15];
permuteidx =  [permuteidx,permuteidx+1];
disp('-- permuted --');
for ii = 1:length(permuteidx)
    disp([num2str(ii),':',fnames(permuteidx(ii)).name]);
end

%% indices
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
for ii = 1:size(permuteidx,2)
    disp([num2str(ii),':',fnames(permuteidx(ii)).name]);
end
MarkProf = MarkProf(:,permuteidx);
MarkProfAvg = (MarkProf(:,1:8)+MarkProf(:,9:16))/2;

%% ----------------------------------------------------------------------------------------
%% --------------------------------         plots          --------------------------------
%% ----------------------------------------------------------------------------------------
%%
incr = 1;
dnameseries = cell(size(MarkProf,2)-2,1); % -2 for control
modseries = {'H3K9me2','H3K27me3'};
repseries = {'rep1','rep2'};
doxseries = {'nodox','dox5days','irreversible','reactivated'};
for ii = 1:length(modseries)
    for kk = 1:length(doxseries)
        for jj = 1:length(repseries)
            disp([doxseries{kk},' ',modseries{ii},' ',repseries{jj}]);
            dnameseries{incr} = [doxseries{kk},' ',modseries{ii},' ',repseries{jj}];
            incr = incr + 1;
        end
    end
end


%% 
xax = -75:5:75; % xax for 5kb
 xax = -77:1:77; % xax for 1kb

SizedFig(50,80);
% subplot(3,1,1);
for ii = 1:length(dnameseries)
    subplot(length(dnameseries)/2,2,ii);
    for jj = 1:length(fnames)
        fnamesp = split(fnames(jj).name,'_');
        dname = [fnamesp{2},' ',fnamesp{3},' ',fnamesp{4}];
        if strcmp(dname,dnameseries{ii})
            bar(xax',flipud(MarkProf(:,jj)),1,'EdgeAlpha',0,'FaceAlpha',1);
            title(dname);
            if strcmp(fnamesp{3},'H3K9me2')
                ylim([0,1.2]);
            elseif strcmp(fnamesp{3},'H3K27me3')
                ylim([0,15]);
%                 ylim([0,1.2]);
            end
        end
    end
%     xlim([-32.5,32.5]);
end
% subplot(3,1,2);
% plot(repmat(xax',[1 7]),MarkProfH3K9me(:,8:14),'o-');
% subplot(3,1,3);
% plot(xax,MarkProf(:,11),'o-');