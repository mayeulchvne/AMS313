function [ errlist, selected_mu, PP, Arb_decomp, Lrb_decomp, ...
           Rmat, f_parallel, f_perp] = RB_greedy_estim_stabilized( tol, Nmax, Xi, BB, AA_decomp, LL_decomp ) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RB_offline :
% Phase offline de la methode base reduite avec residu stabilise
%   
% INPUT * tol: tolerance sur l'erreur d'approximation base reduite
%       * Nmax: taille maximale de la base reduite
%       * Xi : espace des parametres (taille n_train x 2)
%       * BB: matrice du produit scalaire (taille NbDof x NbDof)
%       * AA_decomp : cellarray des matrices de l'operateur (Qa matrices taille NbDof x NbDof)
%       * LL_decomp : cellarray des vecteurs second membre (taille NbDof x 1)
%
% OUTPUT - errlist: (N elements) courbe de convergence
%        - selected_mu: taille Nx2, parametres choisis par l'algo Greedy.
%        - PP: matrice de changement de base (taille NbDof x N)
%        - Arb_decomp: cellarray des parties (taille N x N) de
%                      l'operateur reduit
%        - Lrb_decomp: cellarray des parties (taille N x 1) du second membre reduit 
%        - Rmat: matrice N*Qa x N*Qa du residu stabilise
%        - f_parallel:  matrice N*Qa x Ql du residu stabilise
%        - f_perp:  matrice Ql x Ql du residu stabilise
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nb points de discretisation espace des parametres:
nb_points = size(Xi, 1);
% dimension espace element finis:
NbDof = length(LL_decomp{1});
% nombre de termes decomposition de l'operateur
Qa = length(AA_decomp);
% nombre de termes decomposition du second membre
Ql = length(LL_decomp);

% taille de la base reduite
N = 0;
% erreur d'approximation base reduite
maxestim = 1e9;
% liste des erreur (courbe de convergence)
errlist = [];
% matrice de changement de base
PP = [];
% modele reduit
Arb_decomp = cell(1,Qa);
Lrb_decomp = cell(1,Ql);

% parties du residu independentes du parametre
% --------------------------------------------
% base XX-orthonormale du sous espace XX^{-1}A_q p_n, q=1,...,Qa et
% n=1,...,N (taille Nbpt x N*Qa)
Qmat = [];
% matrice de changement de base (taille N*Qa x N*Qa)
Rmat = [];
% partie du RHS qui appartient au sous espace Range(Qmat), taille Qa*N x Ql
f_parallel=[]; 

% Precalcul des representant de Riesz associes au second membre:
% ------------------
LL_hat_decomp = cell(Ql,1);
for q=1:Ql
    LL_hat_decomp{q} = BB \ LL_decomp{q};
end

% initialisation
% --------------
% On commence au milieu de l'espace des parametres
ind_mu_star = ceil( nb_points/2 );
mu_star = Xi(ind_mu_star, :);
selected_mu=[];

% algorithme Greedy
% ------------------
while (N < Nmax)&&(maxestim > tol)
    
    % sauvegarde valeur du parametre 
    % -------------------------------
    selected_mu = [selected_mu; mu_star];
    
    % calcul de la solution FE avec ce parametre
    % -------------------------- 
    % resolution UU(mu_star)
    UU = PARAMETRIC_solve( mu_star, AA_decomp, LL_decomp );
    
    % orthonormalisation et ajout dans la base
    % ----------------------------------------
    new_basis_func = LINALG_orthonormalize(UU,PP, BB);
    % ajouter la nouvelle fonction de base a la base reduite
    PP = [PP new_basis_func];
    % la taille de la base reduite est augmentee
    N = N + 1;
    disp(sprintf('Reduced basis size N = %i',N));
    
    % mises a jour du solveur RB
    % ---------------------------
    [Arb_decomp, Lrb_decomp] = RB_reduced_decomp(AA_decomp, LL_decomp, PP);
    [Qmat, Rmat]= RB_update_Qmat_Rmat(N, new_basis_func, BB, AA_decomp, Qmat, Rmat);
    [f_parallel, f_perp] = RB_update_resparts_stabilized(Qa, BB, Qmat, LL_decomp, LL_hat_decomp,f_parallel);
    %
    % Probleme de maximisation
    % ------------------------
    % trouver le parametre mu_star qui maximise le residu RB
    [ maxestim, mu_star ] = RB_argmax_estim_stabilized(Xi, Arb_decomp, Lrb_decomp, Rmat, f_parallel, f_perp);
    disp(sprintf('Current Max Estimator = %e', maxestim));
    errlist = [errlist; maxestim];
end

end

