% ============================================================
% Construction d'une base reduite par la methode POD 
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
% construction du maillage
% -------------------------
nx = 200;
ny = 40;
mesh = MESH_build_cartesian(nx, ny);

[ DofNodes, AA_ref, LL_ref,...
      MM, DDX, DDY, BB, AA_decomp, LL_decomp] = FE_assemblages(mesh);


  
n_trial = 100; %valeurs possibles : %à changer
name_PP = strcat('PP_pod_',int2str(n_trial),'.mat');
load(name_PP)
name_mu_list = strcat('mu_list_',int2str(n_trial),'.mat');
load(name_mu_list)



err_br = zeros(n_trial,1);

for i=1:n_trial
    mu = mu_list(i,:);
    %solution haute fidélité
    UU_HF = PARAMETRIC_solve( mu, AA_decomp, LL_decomp);

    %solution base réduite
    [ Arb_decomp, Lrb_decomp ] = RB_reduced_decomp(AA_decomp, LL_decomp, PP);
    Xrb = RB_solve(mu, Arb_decomp, Lrb_decomp);
    Urb = PP*Xrb;
    
    % calcul de l'erreur
    diff = UU_HF - Urb;
    err_br(i) = norm(diff);
end

err_br_max = max(err_br);
err_br_moy = mean(err_br);

fprintf('Erreur base réduite maximale = %e\n',err_br_max)
fprintf('Erreur base réduite moyenne = %e\n',err_br_moy)



