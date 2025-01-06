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
dname = {'DNMT3B nodox','DNMT3B dox12days','DNMT3B/A6 nodox','DNMT3B/A6 dox12days',...
    'KRAB nodox rep1','KRAB dox5days rep1','KRAB irreversible rep1',...
    'KRAB nodox rep2','KRAB dox5days rep2','KRAB irreversible rep2',...
    'Abby control','KRAB reactivated rep1','KRAB reactivated rep2',...
    'KRAB dox5days release2days rep1','KRAB dox5days release5days rep1','KRAB dox5days release10days rep1',...
    'KRAB dox5days release2days rep2','KRAB dox5days release5days rep2','KRAB dox5days release10days rep2',...
    'HDAC4 reactivated rep1','HDAC4 reactivated rep2','HDAC4 irreversible rep1','HDAC4 irreversible rep2',...
    'DNMT3B irreversible rep1','DNMT3B irreversible rep2',...
    'HDAC4 nodox rep1','HDAC4 nodox rep2','HDAC4 dox5days rep1','HDAC4 dox5days rep2'};
dir1idx = [1,2,3,4,5,6,7,8,9, 10,11,NaN,NaN,nan(1,16)];
dir2idx = [1,2,3,4,5,6,7,8,9,NaN,10, 11, 12,nan(1,16)];
dir3idx = [nan(1,13),1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
%% import EM-seq results, heatmap met

% -- 1:unmet, 0:met -- %

cpgidx = [31  ,38  ,51  ,72  ,84  ,103 ,105 ,129 ,142 ,153 ,170 ,189 ,192 ,198 ,208 ,213 ,222 ,243 ,256 ,258 ,275 ,291 ,324 ,339 ,346 ,371 ,381 ,398 ,403 ,430 ,439 ,442 ,444 ,448 ,466 ,468 ,477 ,486 ,527 ,530 ,563 ,590 ,603 ,605 ,609 ,612 ,615 ,622 ,626 ,634 ,644 ,647 ,652 ,661 ,665 ,667 ,674 ,681 ,685 ,706 ,729 ,731 ,734 ,737 ,744 ,749 ,758]+1;

dnamet = {};
for ii = 14%1:length(dir1idx)
    tmp1 = [];
    if ~isnan(dir1idx(ii))
        tmp1 = readmatrix([DIR1,dnames{dir1idx(ii)}] ,'FileType', 'text');
        tmp1 = tmp1(2:end,2:end);
        tmp1 = tmp1(:,cpgidx);
    end
    tmp2 = [];
    if ~isnan(dir2idx(ii))
        tmp2 = readmatrix([DIR2,dnames2{dir2idx(ii)}] ,'FileType', 'text');
        tmp2 = tmp2(2:end,2:end);
        tmp2 = tmp2(:,cpgidx);
    end
    tmp3 = [];
    if ~isnan(dir3idx(ii))
        tmp3 = readmatrix([DIR3,dnames3{dir3idx(ii)}] ,'FileType', 'text');
        tmp3 = tmp3(2:end,2:end);
        tmp3 = tmp3(:,cpgidx);
    end
    tmp = cat(1,tmp1,tmp2,tmp3);
    disp([size(tmp1,1),size(tmp2,1),size(tmp3,1),size(tmp,1)]);
    dnamet = [dnamet,tmp];
end
%% mask unread region 
idxnotread = 27:40;

for ii = 1:length(dnamet)
    tmp = dnamet{ii};
    tmp(:,idxnotread) = NaN;
    dnamet{ii} = tmp;
end

%% write as csv

outDIR = '/Users/tfuji/Documents/tempanalysis/EM-seq/matrices/';
for ii = [1,2,5:10,12:length(dnamet)]
    tmp = dnamet{ii};
    fname = ['EMseq_',strrep(dname{ii},' ','_'),'.csv'];
    writematrix(tmp,[outDIR,fname]);
end

%%
%% bedgraphs
%%
%%
upstreamstart = 1;
upstreamend = 55115765;

reporterstart = 55115766;
reporterend = 55120420;

downstreamstart = 55120421;
downstreamend   = 58622270;

%%
ampliconstart = downstreamstart - 1368;

cpgidx = [31  ,38  ,51  ,72  ,84  ,103 ,105 ,129 ,142 ,153 ,170 ,189 ,192 ,198 ,208 ,213 ,222 ,243 ,256 ,258 ,275 ,291 ,324 ,339 ,346 ,371 ,381 ,398 ,403 ,430 ,439 ,442 ,444 ,448 ,466 ,468 ,477 ,486 ,527 ,530 ,563 ,590 ,603 ,605 ,609 ,612 ,615 ,622 ,626 ,634 ,644 ,647 ,652 ,661 ,665 ,667 ,674 ,681 ,685 ,706 ,729 ,731 ,734 ,737 ,744 ,749 ,758]+1;

cpgoistions = ampliconstart - cpgidx;
cpgoistions_rev = fliplr(cpgoistions);
%%
outDIR = '/Users/tfuji/Documents/tempanalysis/EM-seq/bedgraphs/';
for ii = [1,2,5:10,12:length(dnamet)]
    fname = ['EMseq_',strrep(dname{ii},' ','_'),'.bedgraph'];
    fid = fopen([outDIR,fname],'w');
%     fprintf(fid,'%s\t%d\t%d\t%.2f\n','chr19',0,upstreamend,0.0);
%     fprintf(fid,'%s\t%d\t%d\t%.2f\n','chr19',upstreamend,cpgoistions(end)-1,0.0);
    tmp = dnamet{ii};
    tmp(tmp == -1) = NaN;
    avgtmp = 1 - nanmean(tmp,1);
    avgtmp_rev = fliplr(avgtmp);
    for kk = 1:length(cpgoistions_rev)
        if ~isnan(avgtmp_rev(kk))
         fprintf(fid,'%s\t%d\t%d\t%.2f\n','chr19',...
             cpgoistions_rev(kk),cpgoistions_rev(kk),avgtmp_rev(kk));
        end
    end
    fclose(fid);
end


%%
%%
%%
%% import EM-seq results, heatmap met

% -- 1:unmet, 0:met -- %

cpgidx = [31  ,38  ,51  ,72  ,84  ,103 ,105 ,129 ,142 ,153 ,170 ,189 ,192 ,198 ,208 ,213 ,222 ,243 ,256 ,258 ,275 ,291 ,324 ,339 ,346 ,371 ,381 ,398 ,403 ,430 ,439 ,442 ,444 ,448 ,466 ,468 ,477 ,486 ,527 ,530 ,563 ,590 ,603 ,605 ,609 ,612 ,615 ,622 ,626 ,634 ,644 ,647 ,652 ,661 ,665 ,667 ,674 ,681 ,685 ,706 ,729 ,731 ,734 ,737 ,744 ,749 ,758]+1;

dnamet = {};
for ii = 5%1:length(dir1idx)
    tmp1 = [];
    if ~isnan(dir1idx(ii))
        tmp1 = readmatrix([DIR1,dnames{dir1idx(ii)}] ,'FileType', 'text');
        tmp1 = tmp1(2:end,2:end);
%         tmp1 = tmp1(:,cpgidx);
    end
    tmp2 = [];
    if ~isnan(dir2idx(ii))
        tmp2 = readmatrix([DIR2,dnames2{dir2idx(ii)}] ,'FileType', 'text');
        tmp2 = tmp2(2:end,2:end);
%         tmp2 = tmp2(:,cpgidx);
    end
    tmp3 = [];
    if ~isnan(dir3idx(ii))
        tmp3 = readmatrix([DIR3,dnames3{dir3idx(ii)}] ,'FileType', 'text');
        tmp3 = tmp3(2:end,2:end);
%         tmp3 = tmp3(:,cpgidx);
    end
    tmp = cat(1,tmp1,tmp2,tmp3);
    disp([size(tmp1,1),size(tmp2,1),size(tmp3,1),size(tmp,1)]);
    dnamet = [dnamet,tmp];
end