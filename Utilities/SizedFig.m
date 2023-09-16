function h = SizedFig(PercentOfXwindowSize,PercentOfYwindowSize,Monitor)
    if ~exist('Monitor', 'var'); Monitor = 1; end
    monitorprop = get(groot);
    scrsz = monitorprop.MonitorPositions(Monitor,:);
    figx = scrsz(3)*PercentOfXwindowSize/100;
    figy = scrsz(4)*PercentOfYwindowSize/100;
    xoff = scrsz(3)-figx + scrsz(1) - 1;
    yoff = scrsz(4)-figy + scrsz(2) - 1;
    h = figure('Position',[xoff yoff figx figy]);
    set(h,'PaperUnits','points','PaperSize',[figx, figy]);
end
