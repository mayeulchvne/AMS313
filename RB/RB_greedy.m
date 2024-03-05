function [ errlist, selected_mu, PP, Arb_decomp, Lrb_decomp, ...
           Respart_ll, Respart_la, Respart_aa] = RB_greedy( tol, Nmax, Xi, BB, AA_decomp, LL_decomp ) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RB_offline :
% Phase offline de la methode base reduite
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
%        - Respart_ll: cellarray des scalaires des interations L-L
%        - Respart_la:  cellarray des vecteurs (taille N x 1) des interations L-A
%        - Respart_aa: cellarray des matrices (taille N x N) des interations A-A
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
maxres = 1e9;
% liste des erreur (courbe de convergence)
errlist = [];
% matrice de changement de base
PP = [];
% modele reduit
Arb_decomp = cell(Qa,1);
Lrb_decomp = cell(Ql,1);

% parties du residu independentes du parametre
% --------------------------------------------
% interaction rhs-rhs (termes scalaires)
Respart_ll = cell(Ql,Ql);
% interaction rhs-lhs (des vecteurs taille N x 1)
Respart_la = cell(Ql,Qa);
% interaction rhs-rhs (des matrices taille N x N)
Respart_aa = cell(Qa,Qa);

% Precalcul des representant de Riesz associes au second membre:
% ------------------
Respart_ll = RB_compute_respart_ll( LL_decomp, BB );

% initialisation
% --------------
% On commence au milieu de l'espace des parametres
ind_mu_star = ceil( nb_points/2 );
mu_star = Xi(ind_mu_star, :);
selected_mu=[];

% algorithme Greedy
% ------------------
while (N < Nmax)&&(maxres > tol)
    
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
    [Respart_la, Respart_aa] = RB_update_resparts( PP,BB, AA_decomp,...
                                        LL_decomp, Respart_la, Respart_aa);
    %
    % Probleme de maximisation
    % ------------------------
    % trouver le parametre mu_star qui maximise le residu RB
    [maxres, mu_star] = RB_argmax_resnorm(Xi, Arb_decomp, Lrb_decomp, ...
                                       Respart_ll, Respart_la, Respart_aa);
    disp(sprintf('Current Max Residual = %e', maxres));
    errlist = [errlist; maxres];
end

end

