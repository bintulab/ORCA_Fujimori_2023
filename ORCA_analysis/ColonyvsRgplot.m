%% calculate average angle and Rg

davg = [];
% idxrange = 10:13;
idxrange = 10:13;

davg.angle = [];
davg.Rg = [];
davg.RgRatio = [];
davg.Dist = [];
davg.dname = {};
davg.cfreq_upst = [];
davg.cfreq_dnst = [];

for kk = 1:length(Objs)
    df1 = Objs(kk).df;
    drange = 1:length(df1);
    for ii = 1:length(drange)

        tmp = df1(ii).coordinterp(idxrange,:,:);
        center = round(size(tmp,1)/2);
        upstrm = mean(tmp(1:(center-1),:,:),1) - tmp(center,:,:);
        dnstrm = mean(tmp((center+1):end,:,:),1) - tmp(center,:,:);
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

        tmp = df1(ii).dmatinterp(9:10,12:15,:);
        tmpDist = reshape(mean(mean(tmp,2),1),[],1);
        davg.Dist = [davg.Dist;median(tmpDist)];
        Objs(kk).df(ii).Dist = tmpDist;

        tmp = df1(ii).coordinterp(4:8,:,:);
        tmp = tmp - repmat(mean(tmp,1),[size(tmp,1), 1, 1]);
        tmpRg1 = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
        tmp = df1(ii).coordinterp(8:12,:,:);
        tmp = tmp - repmat(mean(tmp,1),[size(tmp,1), 1, 1]);
        tmpRg2 = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
        tmp = df1(ii).coordinterp(12:16,:,:);
        tmp = tmp - repmat(mean(tmp,1),[size(tmp,1), 1, 1]);
        tmpRg3 = sqrt(reshape(mean(sum(tmp.^2,2),1),[],1));
        davg.RgRatio = [davg.RgRatio;median(tmpRg2./((tmpRg1+tmpRg3)/2))];
        Objs(kk).df(ii).RgRatio = tmpRg2./((tmpRg1+tmpRg3)/2);
        
        tmp = df1(ii).dmatfilt(1:19,1:19,:);
        dnst = zeros(1,9);
        upst = zeros(1,9);
        for jj = 1:9
            tmpupst = reshape(tmp(jj+10,10,:),[],1);
            tmpupst = tmpupst(~isnan(tmpupst));
            dnst(jj) = sum(tmpupst < 100)/length(tmpupst);
            tmpdnst = reshape(tmp(10-jj,10,:),[],1);
            tmpdnst = tmpdnst(~isnan(tmpdnst));
            upst(jj) = sum(tmpdnst < 100)/length(tmpdnst);
        end
        davg.cfreq_dnst = cat(1,davg.cfreq_dnst,dnst);
        davg.cfreq_upst = cat(1,davg.cfreq_upst,upst);
    end
end
% plot
dnameunique = unique(davg.dname);

avgRg = [];
stdRg = [];
avgcfreqU = []; stdcfreqU = [];
avgcfreqD = []; stdcfreqD = [];
for ii = 1:length(dnameunique)
    tmp = davg.Rg(strcmp(davg.dname,dnameunique(ii)));
    avgRg = [avgRg;mean(tmp)];
    stdRg = [stdRg;std(tmp)/sqrt(length(tmp))];
    disp([dnameunique{ii},':',num2str(mean(tmp))]);

    tmp = davg.cfreq_upst(strcmp(davg.dname,dnameunique(ii)),:);
    avgcfreqU = [avgcfreqU;mean(tmp,1)];
    stdcfreqU = [stdcfreqU;std(tmp,[],1)/sqrt(size(tmp,1))];

    tmp = davg.cfreq_dnst(strcmp(davg.dname,dnameunique(ii)),:);
    avgcfreqD = [avgcfreqD;mean(tmp,1)];
    stdcfreqD = [stdcfreqD;std(tmp,[],1)/sqrt(size(tmp,1))];
end

%%
visidx = [1,2,3,4];


% distance
tmpyM = avgRg;%avgAngle;
tmpyS = stdRg;%stdAngle;
ylab = 'median of Rg [nm]';
% cfidx = 9;
% tmpyM = avgcfreqD(:,cfidx);%avgAngle;
% tmpyS = stdcfreqD(:,cfidx);%stdAngle;
% ylab = 'contact freq';
% cfidx = 3;
% tmpyM = avgcfreqU(:,cfidx);%avgAngle;
% tmpyS = stdcfreqU(:,cfidx);%stdAngle;
% ylab = 'contact freq';

rgmin = 0.99*(min(tmpyM(visidx)-tmpyS(visidx))); 
rgmax = 1.01*(max(tmpyM(visidx)+tmpyS(visidx)));

% SizedFig(45,40);
SizedFig(13,30);
% SizedFig(36,30);
disp('---');
for kk = 2
    switch kk
        case 2
            xm = avgcol(visidx); 
            ym = tmpyM(visidx); ys = tmpyS(visidx);
            colmin = -0.5; colmax = 6.5;
            xmin = colmin; xmax = colmax;  ymin = rgmin; ymax = rgmax;
            xlab = 'number of colonies'; ylab = '<median of Rg> [nm]';
            cmap = parula(length(visidx)+2);
            cmap = cmap(2:end-1,:);
    end
    hold on;
    for ii = 1:length(xm)
                plot([xm(ii),xm(ii)],[ym(ii),ym(ii)],'bo','MarkerSize',11,'LineWidth',2,'Color',cmap(ii,:));
        plot([xm(ii)-0,xm(ii)+0]',[ym(ii)-ys(ii),ym(ii)+ys(ii)]','b-','MarkerSize',6,'LineWidth',1.5,...
            'Color',cmap(ii,:),'MarkerEdgeColor',cmap(ii,:),'MarkerFaceColor','w');
%         plot([xm(ii)-(xmax - xmin)/20,xm(ii)+(xmax - xmin)/20],[ym(ii),ym(ii)],'b-','LineWidth',5,'Color',cmap(ii,:));
%         plot([xm(ii)-0,xm(ii)+0]',[ym(ii)-ys(ii),ym(ii)+ys(ii)]','bo','MarkerSize',6,'LineWidth',1,'MarkerEdgeColor',cmap(ii,:),'MarkerFaceColor','w');
        disp(cmap(ii,:));
    end
    xlim([xmin,xmax]);
    ylim([ymin,ymax]);
    if kk == 2
        lm = fitlm(xm,ym);
        round(lm.Rsquared.Ordinary,2)
        intrcpt = lm.Coefficients.Estimate(1);
        slope = lm.Coefficients.Estimate(2);
        plot([xmin,xmax],[xmin,xmax]*slope+intrcpt,'k--','LineWidth',0.5);
    end
    xlabel(xlab); ylabel(ylab);
    set(gca,'FontSize',20,'XColor','k','YColor','k',...
        'XDir','reverse');
    box on;
end



%% contact freq different positions, superimposed
visidx = [1,2,3,4];
visrange = 1:9;%[1,2,5,9]; %
cmapseries = parula(11);
cmapseries =flipud(cmapseries(2:end-1,:));
cmapseries =hot(length(visidx)+10);
cmapseries = cmapseries(1:length(visrange),:);
cmapseries =flipud(cmapseries);%flipud(cmapseries([1,3,7,9],:));

SizedFig(32,32);
disp('---');
ylab = 'contact freq';
R2val = [];
for kk = 1:2
    subplot(1,2,kk);
    cmapseries = fliplr(cmapseries);

    for cfidx = 1:length(visrange)%[1,2,3,4,5,6,7,8,9]%1:9%

            switch kk
                case 1
                    tmpyM = avgcfreqU(:,visrange(cfidx));%avgAngle;
                    tmpyS = stdcfreqU(:,visrange(cfidx));%stdAngle;
                    cmap = cmapseries(cfidx,:);
                    xm = avgcol(visidx); 
                    ym = tmpyM(visidx); ys = tmpyS(visidx);
                    xmin = -0.5; xmax = 6.5;  
                    xlab = {'log col'}; 
                case 2
                    tmpyM = avgcfreqD(:,visrange(cfidx));%avgAngle;
                    tmpyS = stdcfreqD(:,visrange(cfidx));%stdAngle;
                    cmap = cmapseries(cfidx,:);
                    xm = avgcol(visidx); 
                    ym = tmpyM(visidx); ys = tmpyS(visidx);
                    xmin = -0.5; xmax = 6.5;  
                    xlab = {'log col'}; 
            end
            hold on;
            plot(xm,ym,'o-','LineWidth',2,'Color',cmap,'MarkerSize',10,...
                'MarkerFaceColor','w');
            plot([xm-0,xm+0]',[ym-ys,ym+ys]','b-','MarkerSize',6,'LineWidth',1.5,...
        'Color',cmap,'MarkerEdgeColor',cmap,'MarkerFaceColor','w');
            xlim([xmin,xmax]);
            ylim([0.025,0.25]);
            lm = fitlm(xm(2:end),ym(2:end));
            R2val = [R2val;round(lm.Rsquared.Ordinary,2)]
            intrcpt = lm.Coefficients.Estimate(1);
            slope = lm.Coefficients.Estimate(2);
            plot([xmin,xmax],[xmin,xmax]*slope+intrcpt,'--','LineWidth',1,'Color',cmap);
            xlabel(xlab); ylabel(ylab);
            set(gca,'FontSize',20,'XColor','k','YColor','k',...
            'XDir','reverse');
        box on;
    end
end



%% contact freq different positions
visidx = [1,2,3,4];
visrange = 1:9;%[1,2,5,9]; %
cmapseries = parula(11);
cmapseries =flipud(cmapseries(2:end-1,:));
cmapseries =hot(length(visidx)+10);
cmapseries = cmapseries(1:length(visrange),:);
cmapseries =flipud(cmapseries);%flipud(cmapseries([1,3,7,9],:));

SizedFig(32,32);
disp('---');
ylab = 'contact freq';
R2val = [];
for kk = 1:2
    subplot(1,2,kk);
    cmapseries = fliplr(cmapseries);

    for cfidx = 5%:length(visrange)%[1,2,3,4,5,6,7,8,9]%1:9%

            switch kk
                case 1
                    tmpyM = avgcfreqU(:,visrange(cfidx));%avgAngle;
                    tmpyS = stdcfreqU(:,visrange(cfidx));%stdAngle;
                    xm = avgcol(visidx); 
                    ym = tmpyM(visidx); ys = tmpyS(visidx);
                    xmin = -0.5; xmax = 6.5;  
                    xlab = {'log col'}; 
                    cmap = parula(length(visidx)+2);
                    cmap = cmap(2:end-1,:);
                case 2
                    tmpyM = avgcfreqD(:,visrange(cfidx));%avgAngle;
                    tmpyS = stdcfreqD(:,visrange(cfidx));%stdAngle;
                    xm = avgcol(visidx); 
                    ym = tmpyM(visidx); ys = tmpyS(visidx);
                    xmin = -0.5; xmax = 6.5;  
                    xlab = {'log col'}; 
                    cmap = parula(length(visidx)+2);
                    cmap = cmap(2:end-1,:);
            end
            hold on;
            for ii = 1:length(xm)
                plot(xm(ii),ym(ii),'o','LineWidth',2,'Color',cmap(ii,:),'MarkerSize',10,...
                    'MarkerFaceColor','w');
                plot([xm(ii)-0,xm(ii)+0]',[ym(ii)-ys(ii),ym(ii)+ys(ii)]','b-','MarkerSize',6,'LineWidth',1.5,...
            'Color',cmap(ii,:),'MarkerEdgeColor',cmap(ii,:),'MarkerFaceColor','w');
            end
            xlim([xmin,xmax]);
            ylim([0.08,0.23]);
            lm = fitlm(xm(1:end),ym(1:end));
            R2val = [R2val;round(lm.Rsquared.Ordinary,2)]
            intrcpt = lm.Coefficients.Estimate(1);
            slope = lm.Coefficients.Estimate(2);
            plot([xmin,xmax],[xmin,xmax]*slope+intrcpt,'--','LineWidth',1,'Color','k');
            xlabel(xlab); ylabel(ylab);
            set(gca,'FontSize',20,'XColor','k','YColor','k',...
            'XDir','reverse');
        box on;
    end
end


%%
cmap = [0, 0.4470, 0.7410;
       0.8500, 0.3250, 0.0980];

visidx = [1,2,3,4];

% c_freq
disp('---');
for kk = [1:9]
    cfidx = kk;
    tmpyM = avgcfreqU(:,cfidx);%avgAngle;
    tmpyS = stdcfreqU(:,cfidx);%stdAngle;
    xm = avgcol;
    ym = tmpyM(visidx); ys = tmpyS(visidx);
    xmin = -0.1*100; xmax = 1.1*100;  ymin = rgmin; ymax = rgmax;

    lm = fitlm(xm(1:end),ym(1:end));
    RsqU(kk) = round(lm.Rsquared.Ordinary,2);
    intrcpt = lm.Coefficients.Estimate(1);
    slope = lm.Coefficients.Estimate(2);
end
for kk = [1:9]
    cfidx = kk;
    tmpyM = avgcfreqD(:,cfidx);%avgAngle;
    tmpyS = stdcfreqD(:,cfidx);%stdAngle;
    xm = avgcol;
    ym = tmpyM(visidx); ys = tmpyS(visidx);
    xmin = -0.1*100; xmax = 1.1*100;  ymin = rgmin; ymax = rgmax;

    lm = fitlm(xm(1:end),ym(1:end));
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

