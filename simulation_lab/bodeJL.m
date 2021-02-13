% Creating bode plots using the code from listing 1:
% Input variables: G11 = plant, plantTitle = string with title, holdfig
% is an optional parameter which is used to create multiple Bode plots in
% the same figure.

function bodeJL(G11, plantTitle)%, holdfig, figureNr)
   
w = logspace(-4,4,1000)';
Nw = length(w);
[re, im] = nyquist(G11,w);
re = reshape(re, Nw, 1);
im = reshape(im, Nw, 1);
mg11 = 20*log10(sqrt(re.^2+im.^2));
pg11 = 180*phase(re+1i*im)/pi;

% if ~exist('holdfig','var')
%  % holdfig does not exist, so default it to something
%   holdfig = 0;
% end
% 
% if holdfig == 0 % Then create a new figure
%     if ~exist('figureNr','var')
%         % figureNr does not exist, so don't use it.
%         figure();clf;
%     else
%         figure(figureNr);
%     end
% else
%     hold on;
% end
figure();clf;
subplot(2,1,1);
semilogx(w/(2*pi), mg11, 'linewidth', 2); grid;
ylabel('|\cdot|'); 
h = legend(strcat('|', plantTitle, '|')); set(h,'FontSize',18); 
xt = get(gca, 'XTick'); set(gca, 'FontSize', 18)

subplot(2,1,2);
semilogx(w/(2*pi), pg11, 'linewidth', 2); grid;
h = legend(strcat('\phi(',plantTitle,')')); set(h,'FontSize',18); 
xlabel('\omega_[_H_z_]'); ylabel('\phi(\cdot)');
xt = get(gca, 'XTick'); set(gca, 'FontSize', 18)

hold off;
end
