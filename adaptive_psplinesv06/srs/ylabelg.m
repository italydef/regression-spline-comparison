function ylabelg(string)
% YLABEL	Y-axis labels for 2-D and 3-D plots.
% YLABEL('text') adds greek text below the Y-axis 
%   on the current axis.
%   There is a mapping between the standard 
%   keyboard and the symbol keyboard, for example,
%   ylabelg('a') prints an alpha, ylabelg('A')
%   prints a capital alpha.

h = get(gca,'ylabel');
if isempty(h)
	h = text('HorizontalAlignment','center');
	set(gca,'ylabel',h);
end
set(h,'string',string,'FontName','symbol');
