function [dotx,doty] = GENdotPlot(dots,label,ord)

% ul = unique(label);
ul = ord;

dotx = [];doty = [];
for i = 1:length(ul)
    doty = [doty,dots(strcmp(label,ul(i)))];
    dotx = [dotx,i*ones(size(dots(strcmp(label,ul(i)))))];
end