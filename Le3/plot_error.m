function plot_error(figno,data,xhatEKF,xhatEKF_opt,saveplots,filename)

    figure(figno); clf; set(gcf,'WindowStyle','docked');
    
    tiledlayout(2,3,'TileIndexing','columnmajor');
    
    % error signals 
    error.EKF = xhatEKF - data.x;
    error.EKF_opt = xhatEKF_opt - data.x;

    for ii = 1:height(error.EKF)

        nexttile();
        plot( data.t', [xhatEKF(ii,:)', xhatEKF_opt(ii,:)']); hold on
        plot( data.t', data.x(ii,:)','--','Color','k');
        xlabel('$t$ [s]','Interpreter','latex');
        switch ii
            case 1
                ylabel('$x$ [m]','Interpreter','latex')
            case 2
                ylabel('$y$ [m]','Interpreter','latex')
            case 3
                ylabel('$\theta$ [m]','Interpreter','latex')
        end
        box off
        
        nexttile();
        plot( data.t, [error.EKF(ii,:); error.EKF_opt(ii,:)])
        xlabel('$t$ [s]','Interpreter','latex');
        switch ii
            case 1
                ylabel('$\hat x$ [m]','Interpreter','latex')
            case 2
                ylabel('$\hat y$ [m]','Interpreter','latex')
            case 3
                ylabel('$\hat \theta$ [m]','Interpreter','latex')
        end
        box off
    
    end


    % Get latex font in ticks
    h = findall(gcf,'Type','axes'); % An array if you have subplots
    set(h, 'TickLabelInterpreter', 'latex')

    newcolors = ["#1f77b4" "#ff7f0e" "#2ca02c" "#d62728"];
    colororder(newcolors)

    if saveplots
        % saving plots
        textwidth = 14;
        golden_ratio = (1 + sqrt(5)) / 2;
        textheight = textwidth / golden_ratio;
        figsize = [textwidth, textheight];
        
        % Set size and no crop
        set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', figsize);
        set(gcf, 'PaperUnits', 'normalized', 'PaperPosition', [0, 0, 1, 1]);

        savefilename = ['../doc/figures/',filename,'.pdf'];
        print(savefilename,'-dpdf')
    end
end