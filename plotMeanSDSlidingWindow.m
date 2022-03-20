function [fig] = plotMeanSDSlidingWindow(x,y,winLength,NoPoints,maxX, SNReff)
    fig = figure;
    xplot = linspace(min(x),maxX,NoPoints);
    meany = zeros(size(xplot)); SDy = meany;
    for ii = 1:length(xplot)
        Minx = xplot(ii)-winLength/2;
        Maxx = xplot(ii)+winLength/2;
        Idx = (x>=Minx) & (x<Maxx);
        meany(ii) = mean(y(Idx));
        SDy(ii) = std(y(Idx));
    end
    %errorbar(xplot,meany,SDy)
    shadedErrorBar(xplot,meany,SDy,'lineProps','b');
    
    % Quandratic fit
    p = polyfit(xplot,meany,2);
        
    yfit = polyval(p,xplot); hold on;
    plot(xplot,yfit,'r-','linewidth',2);
    text(min(x),max(meany),['polyfit p = [' num2str(p) ']']);
    xlabel kr; ylabel w; title(['w vs. kr, SNReff=', num2str(SNReff)]);
end