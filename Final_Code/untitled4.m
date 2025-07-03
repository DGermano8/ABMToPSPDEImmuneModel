figure;
hold on;
x = linspace(P_A_vals(1), P_A_vals(end), length(P_A_vals));
y = linspace(P_D_vals(1), P_D_vals(end), length(P_D_vals));
[X, Y] = meshgrid(P_A_vals, P_D_vals);

legend_string = [];

cmap = myColour3Gradient(255,yellow_colour,red_colour,purple_colour);


% for kk=1:length(C_sens_vals)
for kk=1:length(C_sens_vals)
    % ModelParams.C_sens = C_sens_vals(kk);
    
    % activation_kk = (PDE_Activation(:,:,kk)-ABM_Activation(:,:,kk))';
    % h = contour(X, Y, activation_kk ,[-0.1,-0.1],'LineWidth',3,'color', cmap(kk,:));

    % max(max(ABM_Activation(:),PDE_Activation(:)))
    
    subplot(1,3,kk)
    hold on;
    contour(X, Y, PDE_Activation(:,:,kk)', ':' ,'LineWidth',2)

    % surf(X, Y, (PDE_Activation(:,:,kk))' )
    % title('PDE')
    % view(2)
    % surf(X, Y, (PDE_Activation-ABM_Activation)' )
    % xlabel('$\alpha_+$','FontSize',18,'Interpreter','latex');
    % ylabel('$\alpha_-$','FontSize',18,'Interpreter','latex');
    colormap(cmap);
    colorbar;
    % legend(legend_string,'Location','north','FontSize',16,'Interpreter','latex');
    % caxis([0 1])

    % subplot(1,2,2)
    % surf(X, Y, (ABM_Activation(:,:,kk))' )
    contour(X, Y, ABM_Activation(:,:,kk)','-' ,'LineWidth',2)

    % view(2)
    % title('ABM')
    % surf(X, Y, (PDE_Activation-ABM_Activation)' )
    % xlabel('$\alpha_+$','FontSize',18,'Interpreter','latex');
    % ylabel('$\alpha_-$','FontSize',18,'Interpreter','latex');
    colorbar;
    title("$\chi = " + num2str(C_sens_vals(kk)) + "$",'FontSize',18,'Interpreter','latex');
    % legend(legend_string,'Location','north','FontSize',16,'Interpreter','latex');
    % caxis([0 1])
    
    % surf(X, Y, log(ABM_Activation(:,:,kk))'-log(PDE_Activation(:,:,kk))' )
    
    legend_string = [legend_string; ["$\chi = " + num2str(C_sens_vals(kk)) + "$"]];
end
% surf(X, Y, (PDE_Activation-ABM_Activation)' )
xlabel('$\alpha_+$','FontSize',18,'Interpreter','latex');
ylabel('$\alpha_-$','FontSize',18,'Interpreter','latex');
% colorbar;
% legend(legend_string,'Location','north','FontSize',16,'Interpreter','latex');
colorbar
hold off;


f = gcf;
f.Color = [1 1 1];
% export_fig 1D_NoTaxis_NoLoss.png -m2.5

