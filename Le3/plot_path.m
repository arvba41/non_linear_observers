%% Plot path 
function plot_path(figno,data,xhatEKF,xhatEKF_opt,saveplots,filename)

figure(figno); clf; set(gcf,'WindowStyle','docked');

tiledlayout(1,3)

%%%% plotting direct estimation
nexttile;
plot( data.x(1,:)', data.x(2,:)', 'LineWidth', 1)
hold on
plot( data.x(1,1)', data.x(2,1)', 'bo', ...
  data.x(1,end)', data.x(2,end)', 'rx');
plot( (data.y(1,:).*cos(data.y(2,:)))', (data.y(1,:).*sin(data.y(2,:)))', 'k')
hold off
axis([0 max(data.x(1,:))*1.1 0 max(data.x(2,:))*1.1]);
axis square
xlabel('$x$ [m]','Interpreter','latex');
xlabel('$y$ [m]','Interpreter','latex');
title('Direct','Interpreter','latex');

%%%% plotting the gussed gains EKF
nexttile;
plot( data.x(1,:)', data.x(2,:)', 'LineWidth', 1)
hold on
plot(xhatEKF(1,:)', xhatEKF(2,:)');
plot( data.x(1,1)', data.x(2,1)', 'bo', ...
  data.x(1,end)', data.x(2,end)', 'rx');
hold off
axis([0 max(data.x(1,:))*1.1 0 max(data.x(2,:))*1.1]);
axis square
xlabel('$x$ [m]','Interpreter','latex');
xlabel('$y$ [m]','Interpreter','latex');
title('EKF','Interpreter','latex');

%%%% plotting the optimizaed gains EKF
nexttile;
plot( data.x(1,:)', data.x(2,:)', 'LineWidth', 1)
hold on
plot(xhatEKF_opt(1,:)', xhatEKF_opt(2,:)');
plot( data.x(1,1)', data.x(2,1)', 'bo', ...
  data.x(1,end)', data.x(2,end)', 'rx');
hold off
axis([0 max(data.x(1,:))*1.1 0 max(data.x(2,:))*1.1]);
axis square
xlabel('$x$ [m]','Interpreter','latex');
xlabel('$y$ [m]','Interpreter','latex');
title('EKF (opti)','Interpreter','latex');

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