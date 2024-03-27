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
NbNodes = size(DofNodes,1);

  
n_trial = 1000; 
N = 8; %valeurs possibles : 4,8,12,16,20,24,28,32
name_PP = strcat('PP_pod_',int2str(N),'.mat');          %changer en PP_greedy!!!!!
load(name_PP)
mu_list=[];
for j=1:n_trial
    mu = [rand*0.9+0.1,rand*0.9+0.1];
    mu_list = [mu_list; mu];
end

[Respart_ll] = RB_compute_respart_ll( LL_decomp, BB );
[Respart_la, Respart_aa] = RB_compute_respart_la_aa(AA_decomp,LL_decomp, PP, BB);

res = zeros(n_trial,1);
alpha = zeros(n_trial,1);

%%%%%%%%%%%%
plan_kappa = true;

if plan_kappa
    err_br = zeros(n_trial,1);
    for i=1:n_trial
        mu = mu_list(i,:);
        alpha(i) = min(mu(1),mu(2))/2;
        
        %calcul solution base réduite
        [ Arb_decomp, Lrb_decomp ] = RB_reduced_decomp(AA_decomp, LL_decomp, PP);
        Xrb = RB_solve(mu, Arb_decomp, Lrb_decomp);
        Urb = PP*Xrb;
        
        %résidu calcul efficace
        res(i) = sqrt(abs(RB_compute_resnorm2(mu, Xrb, Respart_ll, Respart_la, Respart_aa)));
        
        %calcul erreur norme 
        UU_HF = PARAMETRIC_solve( mu, AA_decomp, LL_decomp);
        diff = UU_HF - Urb;
        err_br(i) = sqrt(diff'*BB*diff);
    end
    eff = res./(alpha.*err_br);
    scatter(mu_list(:,1),mu_list(:,2),50,eff,"fill");
end
    
    
    
    
    
    
    

