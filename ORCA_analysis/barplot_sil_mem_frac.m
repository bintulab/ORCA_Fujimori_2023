%% ------------------------------------------------------------------
%% --------------    silenced fraction from RNAFISH    --------------
%% ------------------------------------------------------------------
%% calculate average mRNA and silenced fraction
davgrna = [];
davgrna.spotnum = [];
davgrna.silfrac = [];
davgrna.dname = {};
for kk = 1:2
    switch kk
       case 1
           df1 = df20220802;
       case 2
           df1 = df20220805;
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
%% group by condition
dnameunique = unique(davgrna.dname);
dnameunique = dnameunique([5,3,4,8,6,7,11,9,10,2,1]);
figure;
avgrna = [];
stdrna = [];
for ii = 1:length(dnameunique)
    tmp = davgrna.silfrac(strcmp(davgrna.dname,dnameunique(ii)));
    avgrna = [avgrna;mean(tmp)];
    stdrna = [stdrna;std(tmp)/sqrt(length(tmp))];
    disp([dnameunique{ii},':',num2str(mean(tmp))]);
end
hold on;
bar(avgrna);
plot([1:length(avgrna);1:length(avgrna)],[avgrna-stdrna,avgrna+stdrna]','k-');

for ii = 1:length(dnameunique)
    tmp = davgrna.silfrac(strcmp(davgrna.dname,dnameunique(ii)));
    plot(linspace(ii-0.1,ii+0.1,length(tmp)),tmp,'bo','MarkerFaceColor','white','MarkerEdgeColor','black');
    disp([dnameunique{ii},':',num2str(mean(tmp))]);
end
% ylim([175,210]);

%% bar plot of silenced fraction
dnameplot = dnameunique([1,9,5,6,2,3,11]);
dnameplot = dnameunique([1,5,6,2,3]);
conditionnum = length(dnameplot);
SizedFig(30,30); 
cmap = parula(conditionnum+1);
cmap = [0.6,0.6,0.6;cmap(2:end,:)];

% dnameplot = dnameunique([1,3,7,9]);
% conditionnum = length(dnameplot);
% cmap = repmat([0.6,0.6,0.6;0.3,0.4,1.0],[3 1]);
SizedFig(30,30); 

for i = 1:conditionnum
    hold on;
    dnameplot(i)
    tmp2 = davgrna.silfrac(strcmp(davgrna.dname,dnameplot(i)));
    bar(i,mean(tmp2)*100,'FaceColor',cmap(i,:),'LineWidth',2);
    plot([i-0.2,i+0.2],tmp2*100,'ko','MarkerFaceColor','w','MarkerSize',8,...
        'LineWidth',2);
end
ylim([0,100]);
set(gca,'FontSize',20,'XTick',[]);
% xlabel('after dox removal [days]');
ylabel('silenced fraction [%]');
% legend(string(rec_days),'Location','eastoutside');
box on;

%% ------------------------------------------------------------------
%% ------------------    memory fraction from flow    ---------------
%% ------------------------------------------------------------------
%% check the range of the value
% dat = Output(1:5);
dat = Output;
for ii = 1:length(dat)
    strtmp = split(dat(ii).Name,[" ","_rep","_dox","_rec","_release",".fcs"]);
    dat(ii).cellline = strtmp{2};
    dat(ii).rep = str2double(strtmp{3});
    dat(ii).dox = str2double(strtmp{4});
    dat(ii).rec = str2double(strtmp{5});
    dat(ii).rel = str2double(strtmp{6});
end

%% calculate silenced fraction
idx = 11; % 11=cit
threshold = 15.5;

% idx = 14; % 14=mch
% threshold = 16;

silFrac = [];
for ii = 1:length(dat)
    tmp = dat(ii);
%         ifp = log(tmp(jj).Data(:,17));
%         silFrac(ii,jj) = sum(log(tmp(jj).Data(ifp > 15,idx)) < threshold)/...
%             size(log(tmp(jj).Data(ifp > 15,idx)),1)*100;
    silFrac.frac(ii) = sum(log(tmp.Data(:,idx)) < threshold)/...
        size(log(tmp.Data(:,idx)),1);
    silFrac.rep(ii) = tmp.rep;
    silFrac.rec(ii) = tmp.rec;
    silFrac.rel(ii) = tmp.rel;
    silFrac.dox(ii) = tmp.dox;
    silFrac.cellline{ii} = dat(ii).cellline;
end
%% bar plot of memory fraction by each condition
rec_days = [   0,   1,   5,       5,       5,       1,         1];
cellline = ["A6","A6","A6","A6-147","A6-150","A6-147","A6-HDAC4"];
doxlevel = [   0,1000,1000,    1000,    1000,    1000,      1000];
% cellline = ["B9","B9","B9","B9-147","B9-150","B9-147"];
visorder = [1,6,4,2,3];
conditionnum = length(visorder);
cmap = parula(conditionnum+1);
cmap = [0.6,0.6,0.6;cmap(2:end,:)];

% rec_days = [   0,   5,       0,       5,         0,         1];
% cellline = ["A6","A6","A6-150","A6-150","A6-HDAC4","A6-HDAC4"];
% doxlevel = [   0,1000,       0,    1000,         0,      1000];
% % cellline = ["B9","B9","B9","B9-147","B9-150","B9-147"];
% conditionnum = length(rec_days);
% visorder = [1,2,3,4,5,6];
% cmap = repmat([0.6,0.6,0.6;0.3,0.4,1.0],[3 1]);

SizedFig(30,30); 
for i = 1:conditionnum
    kk = visorder(i);
    hold on;
    plot_criteria = strcmp(silFrac.cellline,cellline(kk)) &...
        (silFrac.dox == doxlevel(kk)) & (silFrac.rec == rec_days(kk));
    tmp1 = silFrac.rel(plot_criteria);
    tmp2 = silFrac.frac(plot_criteria);
    tmp3 = silFrac.rep(plot_criteria);
    tmp2 = tmp2(tmp1 == 16);
%     tmp3 = tmp3(tmp1 == 16);
    bar(i,mean(tmp2)*100,'FaceColor',cmap(i,:),'LineWidth',2);
    plot([i-0.3,i,i+0.3],tmp2(1:3)*100,'ko','MarkerFaceColor','w','MarkerSize',8,...
        'LineWidth',2);
end
ylim([0,100]);
set(gca,'FontSize',20,'XTick',[]);
% xlabel('after dox removal [days]');
ylabel('memory fraction [%]');
% legend(string(rec_days),'Location','eastoutside');
box on;
