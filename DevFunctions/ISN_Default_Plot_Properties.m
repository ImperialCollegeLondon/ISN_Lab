function [PlotLineWidth,PlotFontSize,PlotAxesLineWidth,PlotMarkerSize,PlotPrintResolution] = ISN_Default_Plot_Properties
%ISN_Default_Plot_Properties Sets the default plot properties used in all
%ISN Lab analysis scripts
%   Changes plot LineWidth and MarkerSize and creates plot variables that
%   are called each time a plot is produced. If the plot looks incorrect
%   remember to use set(0,'defaultLineLineWidth',*newvalue*) etc.
PlotLineWidth = 2.5;
PlotFontSize = 22;
PlotAxesLineWidth = 2;
PlotMarkerSize = 10;
set(0,'defaultLineLineWidth',PlotLineWidth);   % set the default line width to lw
set(0,'defaultLineMarkerSize',PlotMarkerSize); % set the default line marker size to msz
PlotPrintResolution = '-r300';
end

