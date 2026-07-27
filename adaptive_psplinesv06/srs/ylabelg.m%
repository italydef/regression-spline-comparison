function xlabelg(string)
%XLABEL	X-axis labels for 2-D and 3-D plots.
% XLABEL('text') adds greek text below the X-axis 
%   on the current axis.
%   There is a mapping between the standard 
%   keyboard and the symbol keyboard, for example,
%   xlabelg('a') prints an alpha, xlabelg('A')
%   prints a capital alpha.

h = get(gca,'xlabel');
if isempty(h)
	h = text('HorizontalAlignment','center');
	set(gca,'xlabel',h);
end
set(h,'string',string,'FontName','symbol');
