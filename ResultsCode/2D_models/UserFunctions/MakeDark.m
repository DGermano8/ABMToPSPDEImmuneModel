function MakeDark()
    set(gcf,'Color',1/255*[56 56 56])
    
    h=findobj(gcf,'type','axes');
    for k=1:size(h,1)
        % f=get(h(k),'children');
        % xi{k}=get(f,'xdata');
        % yi{k}=get(f,'ydata');
        set(h(k),'Color',1/255*[56 56 56])
        set(h(k),'XColor',1/255*[217 217 217])
        set(h(k),'YColor',1/255*[217 217 217])
        set(h(k),'ZColor',1/255*[217 217 217])
        set(h(k).Title,'Color',1/255*[217 217 217])
    end
    
    colorbars = findall(gcf, 'Type', 'ColorBar');
    for k = 1:length(colorbars)
        colorbars(k).Color = 1/255*[56 56 56]; % Change to red or any desired color
        % colorbars(k).TextColor = 1/255*[217 217 217];
        colorbars(k).EdgeColor = 1/255*[217 217 217];
    end
    
    legends = findall(gcf, 'Type', 'legend');
    for k = 1:length(legends)
        legends(k).Color = 1/255*[56 56 56]; % Change to red or any desired color
        legends(k).TextColor = 1/255*[217 217 217];
        legends(k).EdgeColor = 1/255*[217 217 217];
    end


end