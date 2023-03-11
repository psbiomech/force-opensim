function generateOutputFigures(paselected,paquantl,final)

%GENERATEOUTPUTFIGURES Generate output figures - original
%   Prasanna Sritharan, Mario Andrés Muñoz, August 2020
%
% Based on original scripts written by Mario Andrés Muñoz, 2013


% user settings
addpath('..');
user = getUserScriptSettings();
outpath = user.FEATPATH;


fprintf('Generating output figures.\n');
fprintf('------------------------------------------------\n');


% plot quantiles and PC correlations
for f=1:length(final.labels)

    % get PC info from label
    pcstr = '([^_]+)_(\w+_*\w*)_PC(\d+)';
    toks = regexp(final.labels{f},pcstr,'tokens');
    
    % get PC index in quantile struct, and info for plot labels
    switch toks{1}{1}
        case 'angle'
            analysis = 'ik';
            ytext1 = 'Angle (deg)';
        case 'moment'
            analysis = 'id';
            ytext1 = 'Moment (%BW*HT)';
        case 'force'
            analysis = 'so';
            ytext1 = 'Force (BW)';
    end
    
    fprintf('Generating figure: %s...\n',final.labels{f});
    
    figure(f);
    set(gcf,'Position',[100 100 900 400],'PaperSize',[29.7 21.0]);
    sgtitle(strrep(final.labels{f},'_','-'));
    
    % ********************
    % Quantiles
    
    % data
    data.top.mean = paquantl.(analysis).(toks{1}{2}).top.mean(:,str2double(toks{1}{3}));
    data.top.std = paquantl.(analysis).(toks{1}{2}).top.std(:,str2double(toks{1}{3}));
    data.top.upper = data.top.mean + data.top.std;
    data.top.lower = data.top.mean - data.top.std;
    data.bottom.mean = paquantl.(analysis).(toks{1}{2}).bottom.mean(:,str2double(toks{1}{3}));
    data.bottom.std = paquantl.(analysis).(toks{1}{2}).bottom.std(:,str2double(toks{1}{3}));
    data.bottom.upper = data.bottom.mean + data.bottom.std;
    data.bottom.lower = data.bottom.mean - data.bottom.std;
    
    % plot
    subplot(1,2,1);
    hold on;
    qtl = {'top','bottom'};
    clrs = {'b','r'};
    lsty = {'','--'};
    for d=1:2
        [ha,hbot,htop] = shadedplot(0:100,data.(qtl{d}).lower',data.(qtl{d}).upper',clrs{d});
        alpha(ha(2),0.3);
        set(ha,'HandleVisibility','off');   % prevent legend
        delete(hbot);   % remove lower and upper bound lines automatically created by shaded plot
        delete(htop);
        hold on;    % shadedplot sets hold off
        plot(0:100,data.(qtl{d}).mean,[clrs{d} lsty{d}]);        
    end
    hold off;
    legend('25%Q','75%Q');
    box on;
    xlim([0 100]);
    xlabel('% landing phase');
    ylim('auto');   % temporary
    ylabel(ytext1);
    
    
    % ********************
    % Correlations
    
    % data
    coeffs = paselected.(analysis).(toks{1}{2}).coeff(:,str2double(toks{1}{3}));
    normpc2 = paquantl.(analysis).(toks{1}{2}).normpc2(:,str2double(toks{1}{3}));
    
    % plot
    subplot(1,2,2);
    yyaxis left;
    plot(0:100,coeffs,'k');
    ylabel('PC coefficient');
    set(gca,'YColor','k');
    yyaxis right;
    plot(0:100,normpc2,'k--');
    ylabel('Explained variance');
    set(gca,'YColor','k');
    xlim([0 100]);
    xlabel('% landing phase');
    legend('PC coeff','Explained var');
    
    
    
    % output folder
    figdir = fullfile(outpath,'figures');
    if ~exist(figdir,'dir'), mkdir(figdir); end
       
    % save figures and close
    saveas(gcf,fullfile(figdir,[final.labels{f} '.fig']));
    saveas(gcf,fullfile(figdir,[final.labels{f} '.png']));
    saveas(gcf,fullfile(figdir,[final.labels{f} '.pdf']));
    close(gcf);
    
end
    

fprintf('------------------------------------------------\n');

end
   



