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

  
n_trial = 100; 
N = 8; %valeurs possibles : 4,8,12,16,20,24,28,32
name_PP = strcat('PP_pod_',int2str(N),'.mat');
load(name_PP)
mu_list=[];
for j=1:n_trial
    mu = [rand*0.9+0.1,rand*0.9+0.1];
    mu_list = [mu_list; mu];
end

[Respart_ll] = RB_compute_respart_ll( LL_decomp, BB );
[Respart_la, Respart_aa] = RB_compute_respart_la_aa(AA_decomp,LL_decomp, PP, BB);

err_res = zeros(n_trial,1);

validation_residu = false;

if validation_residu
    for i=1:n_trial
        mu = mu_list(i,:);
        %résidu calcul direct
        theta_a = PARAMETRIC_thetaA(mu);
        theta_l = PARAMETRIC_thetaL(mu);

        Qa = length(theta_a);
        Ql = length(theta_l);
    
        AA = zeros(NbNodes,NbNodes);
        for i=1:Qa
            AA = AA + theta_a(i)*AA_decomp{i};
        end
    
        LL = zeros(NbNodes,1);
        for i=1:Ql
            LL = LL + theta_l(i)*LL_decomp{i};
        end
    
        [ Arb_decomp, Lrb_decomp ] = RB_reduced_decomp(AA_decomp, LL_decomp, PP);
        Xrb = RB_solve(mu, Arb_decomp, Lrb_decomp);
        Urb = PP*Xrb;
    
        APx = AA*PP*Xrb-LL;
        y = BB \ APx;
        res_dir = APx'*y;

        %résidu calcul efficace
        res_eff = RB_compute_resnorm2(mu, Xrb, Respart_ll, Respart_la, Respart_aa);
    
    % calcul de l'erreur
    err_res(i) = abs(res_dir-res_eff);
    end

    err_res_max = max(err_res);
    err_res_moy = mean(err_res);
    
    fprintf('Erreur résidu maximale = %e\n',err_res_max)
    fprintf('Erreur résidu moyenne = %e\n',err_res_moy)
end

%%%%%%%%%%%%
convergence_methode = false;
plan_kappa = false;

if convergence_methode
    NN = [];
    E_max = [];
    E_moy = [];
    for i=1:8
        N = 4*i;
        NN = [NN;N];
        name_PP = strcat('PP_pod_',int2str(N),'.mat');
        load(name_PP)
        err_res = zeros(n_trial,1);
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
            err_res(i) = norm(diff);
        end
        E_max = [E_max; max(err_res)];
        E_moy = [E_moy; mean(err_res)];
    end
    semilogy(NN,E_moy,"-o",NN,E_max,"-o");
elseif plan_kappa  
    for i=1=1:n_trial
        
    end
    F = scatteredInterpolant(mu_list(:,1),mu_list(:,2),err_res);
    X = linspace(0.1,1,100);
    Y = linspace(0.1,1,100);
    Z = abs(real(F({X,Y})));
    h = gca;
    surf(X,Y,Z);
    set(h,'zscale','log');
    %shading interp
end
    
    
    
    
    
    
    
    
    
    

