function [] = saveFigureEps(figureName, fontSize)
%SAVEFIGURE Saves the current figure

% By putting this in a separate function it is easy to change the location
% or disable saving for all figures

    hold on;
    if nargin > 1
        set(findobj(gcf,'type','axes'),'FontSize', fontSize);
    else
        set(findobj(gcf,'type','axes'),'FontSize', 18); % default fontSize
    end
    
    imageLocation = '/Users/bas/Dropbox/Apps/Overleaf/Master Thesis/images';
    
    %saveas(gcf,sprintf('%s/%s', imageLocation, figureName), 'epsc');
    
    % texpad not able to live typeset eps, so try pdf 
    
    % Script written by E Akbas (c) Aug 2010 used to remove margins from
    % PDF
    
    saveTightFigure(gcf, sprintf('%s/%s.pdf', imageLocation, figureName));

end

