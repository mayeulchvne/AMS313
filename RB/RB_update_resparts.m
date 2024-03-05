function [Respart_la, Respart_aa] = RB_update_resparts( PP,BB, AA_decomp ,LL_decomp, Respart_la, Respart_aa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RB_update_resparts :
% Mise a jour des parties independentes du parametre du residu
%          
% INPUT * PP: matrice changmeent de base mise a jour (taille Nbpt x N)
%       * BB: matrice produit scalaire (taille Nbpt x Nbpt)
%       * AA_decomp: cellarray des Qa matrices de l'operteur (taille Nbpt x Nbpt)
%       * LL_decomp: cellarray des Ql vecteurs du second membre (taille Nbpt x 1)
%       * Respart_la: cellarray Ql x Qa des vecteurs (taille (N-1) x 1) des interations L-A
%       * Respart_aa: cellarray Qa x Qa des matrices (taille (N-1) x (N-1)) des interations A-A
%
% OUTPUT - Respart_la: cellarray Ql x Qa des vecteurs (taille N x 1) des interations L-A
%        - Respart_aa: cellarray Qa x Qa des matrices (taille N x N) des interations A-A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% taille de la base reduite mise a jour:
% ---------------------------------------
N = size(PP,2);
% derniere fonction de base ajoutee
new_basis_func = PP(:,N);

% Calcul des representants de Riesz de Aq*new_basis_func
% ------------------------------------------------------
Qa = length(AA_decomp);

AqZ_hat = cell(1,Qa);
for q=1:Qa
    AqZ = AA_decomp{q}*new_basis_func;
    AqZ_hat{q} = BB \ AqZ;
end

% Interactions RHS-LHS
% --------------------
Ql = length(LL_decomp);

for q=1:Ql
    for k=1:Qa
        % ajout de la Nieme composante du vecteur P'*Ak'*X^{-1}*Lq
        Respart_la{q,k} = [Respart_la{q,k}; LL_decomp{q}'*AqZ_hat{k}];
    end
end

% Interactions RHS-RHS
% --------------------

for q=1:Qa
    for k=q:Qa
        % calcul et ajout de la Nieme ligne et colonne de la matrice P'*Aq'*X^{-1}*A2*P
        Nth_col = PP'*AA_decomp{q}'*AqZ_hat{k};
        Respart_aa{q,k} = [Respart_aa{q,k} Nth_col(1:(N-1))];
        Nth_row = (PP'*AA_decomp{k}'*AqZ_hat{q})';
        Respart_aa{q,k} = [Respart_aa{q,k} ; Nth_row];
        % remplissage du reste par symetrie
        Respart_aa{k,q} = Respart_aa{q,k}';
    end
end

end

