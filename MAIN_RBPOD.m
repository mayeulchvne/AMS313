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
% definitions 
% -------------------------
% bornes du domaine parametrique
kappa1_min = 0.1;
kappa1_max = 1;
kappa2_min = 0.1;
kappa2_max = 1;
type_plan_experience = 'cartesien'; % 'cartesien' ou 'random'
sauvegarde_base_reduite = false;    % true pour sauvegarder la base reduite
                                    % PP dans un fichier "PP.mat"

% -------------------------
% construction du maillage
% -------------------------
nx = 200;
ny = 40;
mesh = MESH_build_cartesian(nx, ny);

% assemblages des matrices EF
% ---------------------------
[ DofNodes, AA_ref, LL_ref,...
      MM, DDX, DDY, BB, AA_decomp, LL_decomp] = FE_assemblages(mesh);
  
% nombre de degres de liberte
NbDof = size(DofNodes,1);

disp('--------------------');
disp(' Phase Exploratoire ');
disp('--------------------');
tic;

% definition du plan d'experience 
% --------------------------------
% on construit tous les couples (kappa1, kappa2)
% pour lesquels le probleme EF sera resolu (n_train couples en tout)
mu_list = []; % (taille n_train x 2)
if (strcmp(type_plan_experience, 'cartesien'))
    error('type_plan_experience cartesien not yet implemented')
elseif (strcmp(type_plan_experience,'random'))
    error('type_plan_experience random not yet implemented')
else
    error('type_plan_experience pas bien defini')
end

% resolution EF pour tout les couples (kappa1, kappa2) 
% -------------------------------------------------
n_train = size(mu_list,1);
disp(sprintf('Nbre de resolutions HF = %i', n_train));
AllUUs = zeros(NbDof,n_train);
for j=1:n_train
     % valeur du parametre
    mu_j = mu_list(j,:); 
    % calcul de la solution HF
    UU = PARAMETRIC_solve( mu_j, AA_decomp, LL_decomp );
    % on sauvegarde la solution dans la colonne j de AllUUs
    AllUUs(:,j) = UU;
end

elapsed = toc;
disp(sprintf('Phase Exploration elapsed time = %f s', elapsed));

disp('-------------------');
disp(' Phase Compression ');
disp('-------------------');
tic;


elapsed = toc;
disp(sprintf('Phase Compression elapsed time = %f s', elapsed));

% Courbes de decroissance des valeurs propres
% ------------------------------------------------


% Troncature et definition d'une base reduite
% -------------------------------------------
PP = []; % matrice (taille NbDof x N) representant une base reduite de taille N


% Sauvergarde de la base reduite construite
% -----------------------------------------
if(sauvegarde_base_reduite)
    disp('Sauvegarde de la base reduite dans le fichier PP_pod.mat');
    save('PP_pod.mat', 'PP');
end