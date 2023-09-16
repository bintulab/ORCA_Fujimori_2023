function [h,cbar] = imagetriu(input,dmin,dmax,cmap,bcbar,xticklabel,equalaxis)
    ymax = size(input,1);
    yax = 1:ymax;

    if nargin < 2; dmin = min(input(:)); end
    if nargin < 3; dmax = max(input(:)); end
    if nargin < 4; cmap = flipud(parula(256)); end
    if nargin < 5; bcbar = 1; end
    if nargin < 6; xticklabel = yax; end
    if nargin < 7; equalaxis = 1; end

    theta = 45/90*1/2*pi;
    rotR = [cos(theta),-sin(theta);sin(theta),cos(theta)];
    colidx = round(LinTrans(input,dmin,dmax,1,size(cmap,1)));

    hold on;
    for yidx = yax
        for xidx = yidx:ymax
            shift = rotR*[xidx-1;-(yidx-1)]/sqrt(2);
            if xidx == yidx
                xbox = [-1/2 0.0 1/2 0.0];
                ybox = [ 0.0 0.0 0.0 1/2];
                fill(xbox+1+shift(1),ybox+shift(2),cmap(colidx(yidx,xidx),:),'EdgeColor','none');
            else
                xbox = [-1/2 0.0 1/2 0.0];
                ybox = [ 0.0 -1/2 0.0 1/2];
                fill(xbox+1+shift(1),ybox+shift(2),cmap(colidx(yidx,xidx),:),'EdgeColor','none');
            end
        end
    end
    colormap(cmap);
    cbar = colorbar;
    cbar.LineWidth = 1.5;
    cbar.Ticks = cbar.Ticks(1:2:end);
    cbar.TickLabels = LinTrans(cbar.Ticks,0.0,1.0,dmin,dmax);
    if ~bcbar
        cbar.Visible = 'off';
    end
    toppoint = rotR*[ymax-0.5;0.5]/sqrt(2);
    bboxx = [ymax/2 (ymax-0.5) toppoint(1) -1/2 ymax/2]+1;
    bboxy = [0 0 toppoint(2) 0 0];
    plot(bboxx,bboxy,'k-','LineWidth',1.5);
    if equalaxis
        axis image;
    end
    h = gca;
    h.YAxis.Visible = 'off';
    h.Box = 'off';
    set(gca, 'XTick', yax, 'XTickLabels', xticklabel);
end