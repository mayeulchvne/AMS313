function [f_parallel, f_perp] = RB_update_resparts_stabilized(Qa, BB, Qmat, LL_decomp, LL_hat_decomp, f_parallel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RB_update_Qmat_Rmat :
% Mise a jour de la base orthonormale ou vit le residu
%          
% INPUT * Qa: nombre de termes dans la decomposition du LHS
%       * BB: matrice produit scalaire (taille NbDof x NbDof)
%       * Qmat: taille NbDof x Qa*N (base orthonormale ou vit le residu)
%       * LL_decomp: cellarray des Ql vecteur du RHS (taille NbDof x 1) 
%       * LL_hat_decomp: cellarray des Ql vecteur des representants de Riesz du RHS (taille NbDof x 1) 
%       * f_parallel: a mettre a jour, taille Qa*(N-1) x Ql
%
% OUTPUT - f_parallel: mise a jour, taille Qa*N x Ql
%        - f_perp: taille Ql x Ql
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mise a jour f_parallel (rappel; f_parallel := Qmat'*LL_decomp)
% ---------------------------
Ql = length(LL_decomp);

% on va ajouter Qa nouvelles lignes a f_parallel
f_parallel = [f_parallel; zeros(Qa,Ql)];

for k=1:Ql
    % les Qa nouveaux elements de la colonne k de f_parallel
    new_rows= Qmat(:, (end-Qa+1):end)'*LL_decomp{k};
    % ajout a f_parallel
    f_parallel((end-Qa+1):end, k) = new_rows;
end


% calcul de f_perp (rappel; f_perp(q,k) := y(q)'*BB*y(k) 
%                           avec y(q):= LL_hat_decomp{q} - Qmat*f_parallel(:,q))
% ------------------
Nbpt= size(Qmat, 1);
y=zeros(Nbpt, Ql);

for k=1:Ql
    y(:,k)= LL_hat_decomp{k} - Qmat*f_parallel(:,k);
end

f_perp=zeros(Ql, Ql);
for q=1:Ql
    for k=q:Ql
        f_perp(q,k)= y(:,q)'*BB*y(:,k);
        % symetrie
        f_perp(k,q)= f_perp(q,k)';
    end
end

end

