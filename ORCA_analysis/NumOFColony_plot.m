%%

% DIR = '\\171.64.126.245\BintuLab\Taihei\project_ORCA\mES_Nanog_ORCA\20220906_mESregrow_difday1-3\NumberOFColonyEst\';
DIR = '/Volumes/BintuLab/Taihei/project_ORCA/mES_Nanog_ORCA/20220906_mESregrow_difday1-3/NumberOFColonyEst/';
dname = {'dif 1day','dif 1day','dif 1day','dif 1day','dif 1day','dif 1day',...
    'dif 2days','dif 2days','dif 2days','dif 2days','dif 2days','dif 2days',...
    'dif 0day','dif 0day','dif 0day','dif 0day','dif 0day','dif 0day',...
    'dif 3days','dif 3days','dif 3days','dif 3days','dif 3days','dif 3days'};
err = readmatrix([DIR,'error_outsidethewell.txt']);
dcolony = [];
dcolony.colnum = [];
dcolony.commitnum = [];
dcolony.dname = {};
for ii = 1:24
    res = readmatrix([DIR,'TileScan ',sprintf('%03d',ii+5),'.csv']);
    dcolony.colnum = [dcolony.colnum;(size(res,1) - err(ii))];
    dcolony.commitnum = [dcolony.commitnum;500 - (size(res,1) - err(ii))];
    dcolony.dname = [dcolony.dname;dname{ii}];
end
dcolony.colnum = dcolony.colnum([13:18,1:6,7:12,19:24]);
dcolony.dname = dcolony.dname([13:18,1:6,7:12,19:24]);

%% group by condition

dnameunique = unique(dcolony.dname);
SizedFig(15,20);
semilogy(1,1);
% semilogy([1,1]);
avgcol = [];
stdcol = [];
for ii = 1:length(dnameunique)
    tmp = dcolony.colnum(strcmp(dcolony.dname,dnameunique(ii)));
    avgcol = [avgcol;mean(tmp)];
    stdcol = [stdcol;std(tmp)/sqrt(length(tmp))];
    disp([dnameunique{ii},':',num2str(mean(tmp))]);
end
hold on;
cmap = [0, 0.4470, 0.7410];
bar(avgcol,'FaceColor',cmap);
% plot([1:length(avgcol);1:length(avgcol)],[avgcol-stdcol,avgcol+stdcol]','k-');

for ii = 1:length(dnameunique)
    tmp = dcolony.colnum(strcmp(dcolony.dname,dnameunique(ii)));
    if ii == 4
        tmp([1,2]) = 0.4;
    end
    plot(linspace(ii-0.3,ii+0.3,length(tmp)),tmp,'bo','MarkerFaceColor','white','MarkerEdgeColor','black');
    disp([dnameunique{ii},':',num2str(mean(tmp))]);
end
% ylim([175,210]);
set(gca,'XTick',1:4,'XTickLabel',{'day0','day1','day2','day3'},...
    'FontSize',12,'XColor','k','YColor','k','YTick',[1,10,50,100,200,400]);
box off;
ylabel('number of colonies');
ylim([0.4,500]);

%% group by condition

dnameunique = unique(dcolony.dname);
SizedFig(15,20);
% semilogy([1,1]);
avgcol = [];
stdcol = [];
for ii = 1:length(dnameunique)
    tmp = dcolony.colnum(strcmp(dcolony.dname,dnameunique(ii)));
    avgcol = [avgcol;log(mean(tmp))];
    disp([dnameunique{ii},':',num2str(mean(tmp))]);
end
hold on;
bar(avgcol);
% plot([1:length(avgcol);1:length(avgcol)],[avgcol-stdcol,avgcol+stdcol]','k-');

% ylim([175,210]);
set(gca,'XTick',1:4,'XTickLabel',{'day0','day1','day2','day3'},...
    'FontSize',12);
ylabel('number of colonies');

