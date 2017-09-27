% Script to tile multiple figures into one larger one for easy printing
% NOTE: as of yet, doesn't work well combining figures with multiple
% subplots

% credit for version 1.0: https://www.mathworks.com/matlabcentral/answers/7978-put-figures-fig-in-one-page
% Version: 1.0
% REQUIRES HAVING ALL NEEDED FIGURES OPEN

num_figures = 3;

fh = figure;
for ii = 1:(num_figures+1)
    subplot(2,2,ii)
    P{ii} = get(gca,'pos');
end
clf
F = findobj('type','figure');
for ii = 2:(num_figures+2)
    ax = findobj(F(ii),'type','axes');
    set(ax,'parent',fh,'pos',P{(num_figures+3)-ii})
    close(F(ii))
end