function [ Qmat, Rmat ] = RB_update_Qmat_Rmat(N, new_basis_func, BB, AA_decomp, Qmat, Rmat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RB_update_Qmat_Rmat :
% Mise a jour de la base orthonormale ou vit le residu base reduite
%          
% INPUT * new_basis_func: nouvelle fonction de base (taille NbDof x 1)
%       * BB: matrice produit scalaire (taille NbDof x NbDof)
%       * AA_decomp: cellarray des Qa matrices de l'operteur (taille NbDof x NbDof)
%       * Qmat: taille NbDof x Qa*(N-1)
%       * Rmat: (taille Qa*(N-1) x Qa*(N-1)) 
%
% OUTPUT - Qmat: mise a jour, taille NbDof x Qa*N
%        - Rmat: mise a jour, taille Qa*N x Qa*N
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NbDof= length(new_basis_func);

% Calcul des representants de Riesz de Aq*new_basis_func
% ------------------------------------------------------
Qa = length(AA_decomp);

AqZ_hat = zeros(NbDof,Qa);
for q=1:Qa
    AqZ = AA_decomp{q}*new_basis_func;
    AqZ_hat(:,q) = BB \ AqZ;
end

% on orthonormalise les Qa vecteurs BB^{-1} Aq*new_basis_func
% contre les colonnes de Qmat
% ----------------------------------------------------------
for q=1:Qa
    AqZ_hat_perp= LINALG_orthonormalize( AqZ_hat(:,q) , Qmat , BB );
    % ajout dans Qmat (1 colonne en plus):
    Qmat = [Qmat  AqZ_hat_perp];
    %
    % on recupere les coordonnes de AqZ_hat dans la base Qmat 
    % il y a (Qa*(N-1)+q) coordonnees (le nbre de colonne de Qmat)
    AqZ_hat_coords= Qmat'*BB*AqZ_hat(:,q);
    % ajout dans Rmat:
    % ajout d'une ligne de 0
    Rmat= [Rmat; zeros(1,size(Rmat,1)) ];
    % ajout d'une colonne contenant les coordonnees AqZ_hat_coords
    Rmat= [Rmat  AqZ_hat_coords];
end

end

