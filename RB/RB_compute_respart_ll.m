function [Respart_ll] = RB_compute_respart_ll( LL_decomp, BB )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RB_compute_respart_ll :
% calcule le terme d'interaction rhs-rhs dans le residu
%          
% INPUT * LL_decomp: cellarray des parties second membre (Ql vecteurs de taille NbDof x 1)
%       * BB: matrice du produit scalaire (taille NbDof x NbDof)
%
% OUTPUT - Respart_ll:  cellarray (taille Ql x Ql, contenant des scalaires) 
%            des termes d'interaction rhs-rhs dans le residu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nombre de termes second membres
% -------------------------------
Ql = length(LL_decomp);

% calcul des representants de Risez
% ---------------------------------
LL_hat_decomp = cell(1,Ql);
for q=1:Ql
    LL_hat_decomp{q} = BB \ LL_decomp{q};
end

% calcul des interactions
% -----------------------
Respart_ll = cell(Ql,Ql);
for q=1:Ql
    for k=q:Ql
        Respart_ll{q,k} = LL_decomp{k}'*LL_hat_decomp{q};
        % on complete par symetrie:
        Respart_ll{k,q} = Respart_ll{q,k};
    end
end

end

