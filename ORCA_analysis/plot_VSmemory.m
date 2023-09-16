%% calculate average angle and Rg
davg = [];
% idxrange = [1,2,3,4,7,10,13,16,17,18,19];
% idxrange = 4:16; % [4,5,6,7,8,9,10,11,12,13,14,15,16];
idxrange = 8:12; % [4,5,6,7,8,9,10,11,12,13,14,15,16];
% idxrange = 7:13; % [4,5,6,7,8,9,10,11,12,13,14,15,16];
% idxrange = 6:14; % [4,5,6,7,8,9,10,11,12,13,14,15,16];
% idxrange = [7,8,9,10,11];
% idxrange = [10,11,12,13];

davg.angle = [];
davg.Rg = [];
davg.RgFrac = [];
davg.segN = [];
davg.upstsegN = [];
davg.dnstsegN = [];
davg.d2upst = [];
davg.d2dnst = [];
davg.d2cent = [];
davg.cfreq_upst = [];
davg.cfreq_dnst = [];

davg.RgU = [];
davg.RgD = [];

davg.dname = {};

for kk = 1:2%length(Objs)
    df1 = Objs(kk).df;
    drange = 1:length(df1);
    for ii = 1:length(drange)

        tmp = df1(ii).coordinterp(idxrange,:,:);
        center = round(size(tmp,1)/2);
        upstrm = mean(tmp(1:(center-1),:,:),1) - tmp(center,:,:);
        dnstrm = mean(tmp((center+1):end,:,:),1) - tmp(center,:,:);
%         upstrm = tmp(1,:,:) - tmp(center,:,:);
%         dnstrm = tmp(end,:,:) - tmp(center,:,:);
        unorm = reshape(sqrt(sum(upstrm.^2,2)),[],1);
        dnorm = reshape(sqrt(sum(dnstrm.^2,2)),[],1);
        uddot = reshape(sum(upstrm.*dnstrm,2),[],1);
        tmpAngle = acos(uddot./(unorm.*dnorm))/pi*180;
        davg.angle = [davg.angle;median(tmpAngle)];
        davg.dname = [davg.dname,df1(ii).dname];
        Objs(kk).df(ii).angle = tmpAngle;
        
        tmp = df1(ii).coordinterp(idxrange,:,:);
        tmp = tmp - repmat(mean(tmp,1),[size(tmp,1), 1, 1]);
        tmpRg = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
        davg.Rg = [davg.Rg;median(tmpRg)];
        Objs(kk).df(ii).Rg = tmpRg;
        davg.RgFrac = [davg.RgFrac;sum(tmpRg < 150)/length(tmpRg)];
        
        tmp = df1(ii).coordinterp(idxrange,:,:);
        center = round(size(tmp,1)/2);
        dnstrm = sqrt(sum((tmp(1:(center-1),:,:) - tmp(center,:,:)).^2,2));
        upstrm = sqrt(sum((tmp((center+1):end,:,:) - tmp(center,:,:)).^2,2));
        dnstrm = reshape(sum((dnstrm < 150),1),[],1);
        upstrm = reshape(sum((upstrm < 150),1),[],1);
        davg.upstsegN = [davg.upstsegN;mean(upstrm)];
        davg.dnstsegN = [davg.dnstsegN;mean(dnstrm)];
        segN = sqrt(sum((tmp - tmp(center,:,:)).^2,2));
        segN = reshape(sum((segN < 150),1),[],1);
        davg.segN = [davg.segN;mean(segN)];

%         tmp = df1(ii).coordinterp(idxrange,:,:);
%         dnstrm = sqrt(sum(mean(tmp(1:(center),:,:) - tmp(center,:,:),1).^2,2));
%         upstrm = sqrt(sum(mean(tmp((center):(end),:,:) - tmp(center,:,:),1).^2,2));
%         dnstrm = reshape(dnstrm,[],1);
%         upstrm = reshape(upstrm,[],1);
%         davg.d2upst = [davg.d2upst;median(upstrm)];
%         davg.d2dnst = [davg.d2dnst;median(dnstrm)];

        tmp = df1(ii).coordinterp(idxrange,:,:);
        dnstrm = sqrt(sum((tmp(1:(center-1),:,:) - tmp(center,:,:)).^2,2));
        upstrm = sqrt(sum((tmp((center+1):end,:,:) - tmp(center,:,:)).^2,2));
        dnstrm = reshape(nanmean(dnstrm,1),[],1);
        upstrm = reshape(nanmean(upstrm,1),[],1);
        davg.d2upst = [davg.d2upst;median(upstrm)];
        davg.d2dnst = [davg.d2dnst;median(dnstrm)];

        tmp = df1(ii).coordinterp(idxrange,:,:);
        d2cent = sqrt(sum((mean(tmp,1) - tmp(center,:,:)).^2,2));
        d2cent = reshape(d2cent,[],1);
        davg.d2cent = [davg.d2cent;median(d2cent)];

        tmp = df1(ii).coordinterp(idxrange,:,:);
        tmp = tmp((center):(center+1),:,:); % tmp((center):(center+3),:,:);
        tmp = tmp - repmat(mean(tmp,1),[size(tmp,1), 1, 1]);
        tmpRg = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
        davg.RgU = [davg.RgU;median(tmpRg)];
        
        tmp = df1(ii).coordinterp(idxrange,:,:);
        tmp = tmp((center-1):(center),:,:); % tmp = tmp((center-3):(center),:,:);
        tmp = tmp - repmat(mean(tmp,1),[size(tmp,1), 1, 1]);
        tmpRg = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
        davg.RgD = [davg.RgD;median(tmpRg)];

        tmp = df1(ii).dmatfilt(1:19,1:19,:);
        upst = zeros(1,9);
        dnst = zeros(1,9);
        for jj = 1:9
            tmpupst = reshape(tmp(jj+10,10,:),[],1);
            tmpupst = tmpupst(~isnan(tmpupst));
            upst(jj) = sum(tmpupst < 150)/length(tmpupst);
            tmpdnst = reshape(tmp(10-jj,10,:),[],1);
            tmpdnst = tmpdnst(~isnan(tmpdnst));
            dnst(jj) = sum(tmpdnst < 150)/length(tmpdnst);
        end
        davg.cfreq_upst = cat(1,davg.cfreq_upst,upst);
        davg.cfreq_dnst = cat(1,davg.cfreq_dnst,dnst);
    end
end
%% calculate avg and sem 
% plot
dnameunique = unique(davg.dname);
% dnameunique = dnameunique([5,3,4,8,6,7,11,9,10,13,12,2,1]);
dnameunique = dnameunique([3,1,2,6,4,5,9,7,8]);

avgRg = []; stdRg = [];
avgRgU = []; stdRgU = [];
avgRgD = []; stdRgD = [];
avgRgFrac = []; stdRgFrac = [];
avgsegN = []; stdsegN = [];
avgupstsegN = []; stdupstsegN = [];
avgdnstsegN = []; stddnstsegN = [];
avgd2upst = []; stdd2upst = [];
avgd2dnst = []; stdd2dnst = [];
avgd2cent = []; stdd2cent = [];
avgcfreqU = []; stdcfreqU = [];
avgcfreqD = []; stdcfreqD = [];

for ii = 1:length(dnameunique)
    tmp = davg.Rg(strcmp(davg.dname,dnameunique(ii)));
    avgRg = [avgRg;mean(tmp)];
    stdRg = [stdRg;std(tmp)/sqrt(length(tmp))];
    disp([dnameunique{ii},':',num2str(mean(tmp))]);

    tmp = davg.RgFrac(strcmp(davg.dname,dnameunique(ii)));
    avgRgFrac = [avgRgFrac;mean(tmp)];
    stdRgFrac = [stdRgFrac;std(tmp)/sqrt(length(tmp))];

    tmp = davg.segN(strcmp(davg.dname,dnameunique(ii)));
    avgsegN = [avgsegN;mean(tmp)];
    stdsegN = [stdsegN;std(tmp)/sqrt(length(tmp))];

    tmp = davg.upstsegN(strcmp(davg.dname,dnameunique(ii)));
    avgupstsegN = [avgupstsegN;mean(tmp)];
    stdupstsegN = [stdupstsegN;std(tmp)/sqrt(length(tmp))];

    tmp = davg.dnstsegN(strcmp(davg.dname,dnameunique(ii)));
    avgdnstsegN = [avgdnstsegN;mean(tmp)];
    stddnstsegN = [stddnstsegN;std(tmp)/sqrt(length(tmp))];

    tmp = davg.d2upst(strcmp(davg.dname,dnameunique(ii)));
    avgd2upst = [avgd2upst;mean(tmp)];
    stdd2upst = [stdd2upst;std(tmp)/sqrt(length(tmp))];

    tmp = davg.d2dnst(strcmp(davg.dname,dnameunique(ii)));
    avgd2dnst = [avgd2dnst;mean(tmp)];
    stdd2dnst = [stdd2dnst;std(tmp)/sqrt(length(tmp))];

    tmp = davg.d2cent(strcmp(davg.dname,dnameunique(ii)));
    avgd2cent = [avgd2cent;mean(tmp)];
    stdd2cent = [stdd2cent;std(tmp)/sqrt(length(tmp))];

    tmp = davg.RgU(strcmp(davg.dname,dnameunique(ii)));
    avgRgU = [avgRgU;mean(tmp)];
    stdRgU = [stdRgU;std(tmp)/sqrt(length(tmp))];

    tmp = davg.RgD(strcmp(davg.dname,dnameunique(ii)));
    avgRgD = [avgRgD;mean(tmp)];
    stdRgD = [stdRgD;std(tmp)/sqrt(length(tmp))];

    tmp = davg.cfreq_upst(strcmp(davg.dname,dnameunique(ii)),:);
    avgcfreqU = [avgcfreqU;mean(tmp,1)];
    stdcfreqU = [stdcfreqU;std(tmp,[],1)/sqrt(size(tmp,1))];

    tmp = davg.cfreq_dnst(strcmp(davg.dname,dnameunique(ii)),:);
    avgcfreqD = [avgcfreqD;mean(tmp,1)];
    stdcfreqD = [stdcfreqD;std(tmp,[],1)/sqrt(size(tmp,1))];
end

%% calculate average p-value
davg = [];
idxrange = [1,2,3,4,7,10,13,16,17,18,19];
idxrange = 4:16;

davg.p_angle = [];
davg.p_Rg = [];
davg.dname = {};
for kk = 1:length(Objs)
    df1 = Objs(kk).df;
    drange = 1:length(df1);
    for ii = 1:length(drange)
        ctrlidx = df1(ii).ctrlidx;
        tmp = ranksum(Objs(kk).df(ctrlidx).angle,Objs(kk).df(ii).angle);
        davg.p_angle = [davg.p_angle;tmp];

        tmp = ranksum(Objs(kk).df(ctrlidx).Rg,Objs(kk).df(ii).Rg);
        davg.p_Rg = [davg.p_Rg;tmp];
        
        davg.dname = [davg.dname,df1(ii).dname];
    end
end


%%
visidx = [1,5,6,2,3];
visidx2 = [4,7,8,9];
% visidx = [1,4,7,8,9,5,6,2,3];

% distance
tmpyM = avgRg;%avgAngle;
tmpyS = stdRg;%stdAngle;
ylab = 'median of Rg (nm)';
rgmin = 175; rgmax = 202;
% rgmin = 120; rgmax = 150;
tmpyM = avgd2upst;%avgAngle;
tmpyS = stdd2upst;%stdAngle;
ylab = 'med of avg dist. to reporter';
% rgmin = 210; rgmax = 270;
% tmpyM = avgd2dnst;%avgAngle;
% tmpyS = stdd2dnst;%stdAngle;
% ylab = 'med of avg dist. to reporter';
% rgmin = 210; rgmax = 270;
% tmpyM = avgd2cent;%avgAngle;
% tmpyS = stdd2cent;%stdAngle;
% ylab = 'med of avg dist. to reporter';

% angle
% tmpyM = avgAngle;
% tmpyS = stdAngle;
% ylab = 'median of U/D angle [nm]';

% tmpyM = avgRgFrac*100;%avgAngle;
% tmpyS = stdRgFrac*100;%stdAngle;
% ylab = '% Rg < 150 nm';
% rgmin = 10; rgmax = 35;
% tmpyM = avgRgFrac*100;%avgAngle;
% tmpyS = stdRgFrac*100;%stdAngle;
% ylab = '% Rg < 150 nm';

% tmpyM = avgsegN;%avgAngle;
% tmpyS = stdsegN;%stdAngle;
% ylab = 'avg # of seg. < 150 nm';
tmpyM = avgupstsegN;%avgAngle;
tmpyS = stdupstsegN;%stdAngle;
ylab = 'avg # of seg. < 150 nm';
% rgmin = 0.9; rgmax = 2.1;
% 
tmpyM = avgdnstsegN;%avgAngle;
tmpyS = stddnstsegN;%stdAngle;
ylab = 'avg # of seg. < 150 nm';
% rgmin = 0.9; rgmax = 2.1;

% tmpyM = avgRgD;%avgAngle;
% tmpyS = stdRgD;%stdAngle;
% ylab = '<median of Rg> (nm)';

% c_freq
% tmpyM = avgcfreqD(:,cfidx);%avgAngle;
% tmpyS = stdcfreqD(:,cfidx);%stdAngle;
% ylab = 'contact freq';

rgmin = 0.99*(min(tmpyM([visidx,visidx2])-tmpyS([visidx,visidx2]))); 
rgmax = 1.01*(max(tmpyM([visidx,visidx2])+tmpyS([visidx,visidx2])));

% SizedFig(45,40);
SizedFig(36,30);
disp('---');
for kk = 1:2
    subplot(1,2,kk);
    switch kk
        case 1
            xm = [1,2,3,4]; xs = [0,0,0,0];
            ym = tmpyM(visidx2); ys = tmpyS(visidx2);
            xmin =0.5; xmax = 4.5;  ymin = rgmin; ymax = rgmax;
%             xlab = {'KRAB Y46A no dox','KRAB EEW25AAA no dox','KRAB EEW25AAA dox 1 day','KRAB EEW25AAA dox 5 days'}; 
            xlab = {'controls'}; 
            cmap = repmat([0.6,0.6,0.6],[4,1]);
        case 2
            xm = avgmem(visidx)*100; xs = stdmem(visidx)/sqrt(3)*100;
            ym = tmpyM(visidx); ys = tmpyS(visidx);
            xmin = -0.1*100; xmax = 1.1*100;  ymin = rgmin; ymax = rgmax;
            cmap = parula(length(visidx)+1);
            cmap = [0.6,0.6,0.6;cmap(2:end,:)];
            xlab = 'memory (%)'; 
    end
    hold on;
    for ii = 1:length(xm)
        plot([xm(ii),xm(ii)],[ym(ii),ym(ii)],'bo','MarkerSize',11,'LineWidth',2,'Color',cmap(ii,:));
%         plot([xm(ii)-(xmax - xmin)/20,xm(ii)+(xmax - xmin)/20],[ym(ii),ym(ii)],'b-','LineWidth',5,'Color',cmap(ii,:));
        plot([xm(ii)-0,xm(ii)+0]',[ym(ii)-ys(ii),ym(ii)+ys(ii)]','b-','MarkerSize',6,'LineWidth',1.5,...
            'Color',cmap(ii,:),'MarkerEdgeColor',cmap(ii,:),'MarkerFaceColor','w');
%         plot([xm(ii)-xs(ii),xm(ii)+xs(ii)]',[ym(ii)-0,ym(ii)+0]','b-','LineWidth',2,'Color',cmap(ii,:));
        disp(cmap(ii,:));
    end
    xlim([xmin,xmax]);
    ylim([ymin,ymax]);
    if kk == 2
        lm = fitlm(xm(2:end),ym(2:end));
        round(lm.Rsquared.Ordinary,2)
        intrcpt = lm.Coefficients.Estimate(1);
        slope = lm.Coefficients.Estimate(2);
        plot([xmin,xmax],[xmin,xmax]*slope+intrcpt,'k--','LineWidth',0.5);
    end
    xlabel(xlab); ylabel(ylab);
    set(gca,'FontSize',20,'XColor','k','YColor','k');
    box on;
end

%% contact freq different positions; superimpose multiple
visidx = [1,5,6,2,3];
visidx2 = [4,7,8,9];
% cmapseries = parula(11);
% cmapseries =flipud(cmapseries(2:end-1,:));
% cmapseries = parula(length(visidx)+2);
% cmapseries =hot(length(visidx)+8);
% cmapseries =flipud(cmapseries([1,3,7,9],:));
% cmapseries =fliplr(cmapseries);
cmapseries =hot(11);
cmapseries =fliplr(cmapseries(2:end-1,:));
visrange = [1,2,5,9];%[2,4,6];%1:9;%
SizedFig(32,32);
disp('---');
ylab = 'contact freq';
for cfidx = 1:length(visrange)%1:9%[1,2,5,9]%
    tmpyM = avgcfreqD(:,visrange(cfidx));%avgAngle;
    tmpyS = stdcfreqD(:,visrange(cfidx));%stdAngle;
    cmap = cmapseries(cfidx,:);

    for kk = 1:2
        subplot(1,2,kk);
        switch kk
            case 1
                xm = [1;2;3;4]; xs = [0,0,0,0];
                ym = tmpyM(visidx2); ys = tmpyS(visidx2);
                xmin =0.5; xmax = 4.5;  
    %             xlab = {'KRAB Y46A no dox','KRAB EEW25AAA no dox','KRAB EEW25AAA dox 1 day','KRAB EEW25AAA dox 5 days'}; 
                xlab = {'controls'}; 
            case 2
                xm = avgmem(visidx)*100; xs = stdmem(visidx)/sqrt(3)*100;
                ym = tmpyM(visidx); ys = tmpyS(visidx);
                
                xmin = -0.1*100; xmax = 1.1*100;  
                xlab = {'memory [%]'}; 
        end
        hold on;
        plot(xm,ym,'o-','LineWidth',2,'Color',cmap,'MarkerSize',10,...
            'MarkerFaceColor','w');
        plot([xm-0,xm+0]',[ym-ys,ym+ys]','b-','MarkerSize',6,'LineWidth',1.5,...
    'Color',cmap,'MarkerEdgeColor',cmap,'MarkerFaceColor','w');

        xlim([xmin,xmax]);
%         ylim([0,0.2]);
%         ylim([0.05,0.4]);
        if kk == 2
            lm = fitlm(xm(2:end),ym(2:end));
            round(lm.Rsquared.Ordinary,2)
            intrcpt = lm.Coefficients.Estimate(1);
            slope = lm.Coefficients.Estimate(2);
            plot([xmin,xmax],[xmin,xmax]*slope+intrcpt,'--','LineWidth',1,'Color',cmap);
        end
        xlabel(xlab); ylabel(ylab);
        set(gca,'FontSize',20,'XColor','k','YColor','k');
        box on;
    end
end


%% contact freq vs memory
visidx = [1,5,6,2,3];
visidx2 = [4,7,8,9];
% visidx = [1,4,7,8,9,5,6,2,3];

% c_freq
cfidx = 1;
tmpyM = avgcfreqD(:,cfidx);%avgAngle;
tmpyS = stdcfreqD(:,cfidx);%stdAngle;
ylab = 'contact freq';
% rgmin = 0.12; 
% rgmax = 0.33;

rgmin = 0.99*(min(tmpyM([visidx,visidx2])-tmpyS([visidx,visidx2]))); 
rgmax = 1.01*(max(tmpyM([visidx,visidx2])+tmpyS([visidx,visidx2])));

% SizedFig(45,40);
SizedFig(30,30);
disp('---');
for kk = 1:2
    subplot(1,2,kk);
    switch kk
        case 1
            xm = [1,2,3,4]; xs = [0,0,0,0];
            ym = tmpyM(visidx2); ys = tmpyS(visidx2);
            xmin =0.5; xmax = 4.5;  ymin = rgmin; ymax = rgmax;
%             xlab = {'KRAB Y46A no dox','KRAB EEW25AAA no dox','KRAB EEW25AAA dox 1 day','KRAB EEW25AAA dox 5 days'}; 
            xlab = {'controls'}; 
            cmap = repmat([0.6,0.6,0.6],[4,1]);
        case 2
            xm = avgmem(visidx)*100; xs = stdmem(visidx)/sqrt(3)*100;
            ym = tmpyM(visidx); ys = tmpyS(visidx);
            xmin = -0.1*100; xmax = 1.1*100;  ymin = rgmin; ymax = rgmax;
            cmap = parula(length(visidx)+1);
            cmap = [0.6,0.6,0.6;cmap(2:end,:)];
            xlab = 'memory (%)'; 
    end
    hold on;
    for ii = 1:length(xm)
        plot([xm(ii),xm(ii)],[ym(ii),ym(ii)],'bo','MarkerSize',11,'LineWidth',2,'Color',cmap(ii,:));
%         plot([xm(ii)-(xmax - xmin)/20,xm(ii)+(xmax - xmin)/20],[ym(ii),ym(ii)],'b-','LineWidth',5,'Color',cmap(ii,:));
        plot([xm(ii)-0,xm(ii)+0]',[ym(ii)-ys(ii),ym(ii)+ys(ii)]','b-','MarkerSize',6,'LineWidth',1.5,...
            'Color',cmap(ii,:),'MarkerEdgeColor',cmap(ii,:),'MarkerFaceColor','w');
%         plot([xm(ii)-xs(ii),xm(ii)+xs(ii)]',[ym(ii)-0,ym(ii)+0]','b-','LineWidth',2,'Color',cmap(ii,:));
        disp(cmap(ii,:));
    end
    xlim([xmin,xmax]);
    ylim([ymin,ymax]);
    if kk == 2
        lm = fitlm(xm(2:end),ym(2:end));
        round(lm.Rsquared.Ordinary,2)
        intrcpt = lm.Coefficients.Estimate(1);
        slope = lm.Coefficients.Estimate(2);
        plot([xmin,xmax],[xmin,xmax]*slope+intrcpt,'k--','LineWidth',0.5);
    end
    xlabel(xlab); ylabel(ylab);
    set(gca,'FontSize',20,'XColor','k','YColor','k');
    box on;
end

%%
cmap = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980];

visidx = [1,5,6,2,3];
visidx2 = [4,7,8,9];
% visidx = [1,4,7,8,9,5,6,2,3];

% c_freq
disp('---');
for kk = [1:9]
    cfidx = kk;
    tmpyM = avgcfreqU(:,cfidx);%avgAngle;
    tmpyS = stdcfreqU(:,cfidx);%stdAngle;
    xm = avgmem(visidx)*100; xs = stdmem(visidx)/sqrt(3)*100;
    ym = tmpyM(visidx); ys = tmpyS(visidx);
    xmin = -0.1*100; xmax = 1.1*100;  ymin = rgmin; ymax = rgmax;

    lm = fitlm(xm(2:end),ym(2:end));
    RsqU(kk) = round(lm.Rsquared.Ordinary,2);
    intrcpt = lm.Coefficients.Estimate(1);
    slope = lm.Coefficients.Estimate(2);
end
for kk = [1:9]
    cfidx = kk;
    tmpyM = avgcfreqD(:,cfidx);%avgAngle;
    tmpyS = stdcfreqD(:,cfidx);%stdAngle;
    xm = avgmem(visidx)*100; xs = stdmem(visidx)/sqrt(3)*100;
    ym = tmpyM(visidx); ys = tmpyS(visidx);
    xmin = -0.1*100; xmax = 1.1*100;  ymin = rgmin; ymax = rgmax;

    lm = fitlm(xm(2:end),ym(2:end));
    RsqD(kk) = round(lm.Rsquared.Ordinary,2);
    intrcpt = lm.Coefficients.Estimate(1);
    slope = lm.Coefficients.Estimate(2);
end
SizedFig(40,15);
hold on;
bar([-75,-60,-45,-30:5:-5],fliplr(RsqU),'BarWidth',1);
bar([5:5:30,45,60,75],RsqD,'BarWidth',1);
plot([-80,0],[mean(RsqU) mean(RsqU)],'--','Color',cmap(1,:));
plot([0,80],[mean(RsqD) mean(RsqD)],'--','Color',cmap(2,:));
box on;
set(gca,'FontSize',15,'XColor','k','YColor','k');
%%
figure;
subplot(1,2,1);
hold on;
histogram(RsqU,10,'BinLimits',[0,1]);
ylim([0,4]);
subplot(1,2,2);
histogram(RsqD,10,'BinLimits',[0,1]);
ylim([0,4]);
%% histogram of distance

cmap = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980];
SizedFig(15,20);
hold on;
binnum = 30;
binr = [0,800];

kk = 1;
df1 = Objs(kk).df;
for ii = 1%:2
    switch ii
        case 1
            tmp = df1(1).dmatfilt(1:19,1:19,:);
            upst = tmp(12,10,:);
            upst = reshape(upst(~isnan(upst)),[],1);
            upst = upst(~isnan(upst));
            sum(upst < 150)/length(upst)
        case 2
            tmp = df1(3).dmatfilt(1:19,1:19,:);
            upst = tmp(12,10,:);
            upst = reshape(upst(~isnan(upst)),[],1);
            upst = upst(~isnan(upst));
            sum(upst < 150)/length(upst)
    end
    h1 = histogram(upst,binnum,'BinLimits',binr,'normalization','probability',...
        'EdgeAlpha',0.0,'FaceColor',cmap(ii,:));
    stairs([h1.BinEdges,h1.BinEdges(end)],[h1.Values,h1.Values(end),0],'k-');

    set(gca,'XColor','k','YColor','k');

    xlabel('distance (nm)');
    ylabel('probability');
    box on;
end
plot([150,150],[0,0.12],'k--');

%% cumsum

cmap = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980];
SizedFig(15,20);
hold on;
binnum = 100;
binr = [0,800];

kk = 1;
df1 = Objs(kk).df;
for ii = 1:2
    switch ii
        case 1
            tmp = df1(1).dmatfilt(1:19,1:19,:);
            upst = tmp(12,10,:);
            upst = reshape(upst(~isnan(upst)),[],1);
            upst = upst(~isnan(upst));
            sum(upst < 150)/length(upst)
        case 2
            tmp = df1(3).dmatfilt(1:19,1:19,:);
            upst = tmp(12,10,:);
            upst = reshape(upst(~isnan(upst)),[],1);
            upst = upst(~isnan(upst));
            sum(upst < 150)/length(upst)
    end
    h1 = histogram(upst,binnum,'BinLimits',binr,'normalization','probability',...
        'EdgeAlpha',0.0,'FaceColor',cmap(ii,:));
%     stairs([h1.BinEdges,h1.BinEdges(end)],[h1.Values,h1.Values(end),0],'k-');
    tmp = [h1.Values,h1.Values(end),0];
    stairs([h1.BinEdges,h1.BinEdges(end)],cumsum(tmp),'k-');

    set(gca,'XColor','k','YColor','k');

    xlabel('distance (nm)');
    ylabel('probability');
    box on;
end
plot([150,150],[0,0.08],'k--');






