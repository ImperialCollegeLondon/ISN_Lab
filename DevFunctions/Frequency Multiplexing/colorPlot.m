function colorPlot(x,y,col,num_colour)
%colorPlot plot x-y plot with colour defined by a column vector
%   x : x data
%   y : y data
%   col : colour data (type: double)
%   num_colour : number of possible colour

if ~isrow(x)
    x = x';
end

if ~isrow(y)
    y = y';
end
if ~isrow(col)
    col = col';
end

if range(col) == 0
    cmap = jet(num_colour);
    figure
    plot(x,y,'Color',cmap(col(1)));
else
    colormap jet
    surface([x;x],[y;y],[zeros(1,length(x));zeros(1,length(x))],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
    colorbar;
end
end

