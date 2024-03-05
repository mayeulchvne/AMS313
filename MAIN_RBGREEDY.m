% ============================================================
% Construction d'une base reduite par la methode Greedy 
% pour l'equation de Poisson 2D, avec conditions de
% Neumann et de Dirichlet sur le bord
%  
% avec 2 parametres:
%   * kappa1 (coeff diffusion dans ss-domaine \Omega1)
%   * kappa2 (coeff diffusion dans ss-domaine \Omega2)
%
% Note: le coeff diffusion dans le ss-domaine \Omega0 est fixe
% ============================================================

close all;
clear all;

% -------------------------
% definitions 
% -------------------------
% bornes du domaine parametrique
kappa1_min = 0.1;
kappa1_max = 1;
kappa2_min = 0.1;
kappa2_max = 1;
classical_or_stabilized = 'classical'; % residu 'classical' ou 'stabilized'
type_plan_entrainement = 'cartesien'; % 'cartesien' ou 'random'
sauvegarde_base_reduite = false;    % true pour sauvegarder la base reduite
                                    % PP dans un fichier "PP.mat"
                                    
% construction du maillage
% -------------------------
nx = 200;
ny = 40;
mesh = MESH_build_cartesian(nx, ny);

% assemblages des matrices EF
% ---------------------------
[ DofNodes, AA_ref, LL_ref,...
      MM, DDX, DDY, BB, AA_decomp, LL_decomp] = FE_assemblages(mesh);

% definition du plan d'entrainement
% ---------------------------------
% on construit tous les couples (kappa1, kappa2)
% pour lesquels le probleme RB sera resolu
mu_train = [];
if (strcmp(type_plan_entrainement, 'cartesien'))
    error('type_plan_entrainement cartesien not yet implemented')
elseif (strcmp(type_plan_entrainement,'random'))
    error('type_plan_entrainement random not yet implemented')
else
    error('type_plan_entrainement pas bien defini')
end

% Lancement de l'algo Greedy
% ---------------------------
disp('-------------------')
disp(' Algorithme Greedy')
disp('-------------------')
tol = 1.E-7; 
Nmax= 30;
if (strcmp(classical_or_stabilized, 'classical'))
    [ errlist, selected_mu, PP, Arb_decomp, Lrb_decomp, ...
      Respart_ll, Respart_la, Respart_aa] = RB_greedy( tol, Nmax, mu_train, BB,AA_decomp, LL_decomp); 
elseif (strcmp(classical_or_stabilized, 'stabilized'))
    [ errlist, selected_mu, PP, Arb_decomp, Lrb_decomp, ...
      Rmat, f_parallel, f_perp] = RB_greedy_stabilized(tol, Nmax, mu_train, BB, AA_decomp, LL_decomp );
else
    error('classical_or_stabilized pas bien defini !!')
end

% Affichage courbe de convergence
% -------------------------------
figure_conv = figure;
axes_conv = axes('Parent',figure_conv,'YScale','log','YMinorTick','on',...
    'YMinorGrid','on',...
    'YGrid','on', 'FontSize',16);
box(axes_conv,'on');
hold(axes_conv,'all');
semilogy(errlist , 'Marker','o','LineWidth',2,'LineStyle',':');
xlabel('RB size $N$', 'interpreter', 'latex', 'Fontsize',18);
ylabel('Max Residual Norm', 'interpreter', 'latex', 'Fontsize',18);
title('Greedy convergence curve', 'interpreter', 'latex', 'Fontsize',18);

% affichage des parametres "optimaux" selectionnes par Greedy
% -----------------------------------------------------------
figure_mu_optim = figure;
axes_mu_optim = axes('Parent',figure_mu_optim,'FontSize',16);
box(axes_mu_optim,'on');
hold(axes_mu_optim,'all');
scatter(selected_mu(:,1), selected_mu(:,2), ...
    'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerEdgeColor',[0 0 1],...
    'Marker','square');
xlim([kappa1_min kappa1_max])
ylim([kappa2_min kappa2_max])
xlabel('$\kappa_1$', 'interpreter', 'latex', 'Fontsize',18);
ylabel('$\kappa_2$', 'interpreter', 'latex', 'Fontsize',18);
title('Greedy optimal parameters', 'interpreter', 'latex', 'Fontsize',18);

% Sauvergarde de la base reduite construite
% -----------------------------------------
if(sauvegarde_base_reduite)
    disp('Sauvegarde de la base reduite dans le fichier PP_greedy.mat');
    save('PP_greedy.mat');
end