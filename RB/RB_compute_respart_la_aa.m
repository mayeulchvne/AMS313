function [Respart_la, Respart_aa] = RB_compute_respart_la_aa(AA_decomp,LL_decomp, PP, BB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RB_compute_respart_la_aa :
% calcule les termes d'interaction lhs-rhs et lhs-lhs dans le residu RB
%          
% INPUT * AA_decomp: cellarray des parties de l'operateur (Qa matrice
%                      de taille NbDof x NbDof)
%       * LL_decomp: cellarray des parties second membre (Ql vecteurs de 
%                      taille NbDof x 1)
%       * PP: base reduite  (taille NbDof x N)
%       * BB: matrice du produit scalaire (taille NbDof x NbDof)
%
% OUTPUT - Respart_la:  cellarray (taille Ql x Qa, contenant des vecteurs 
%            de taille N) des termes d'interaction rhs-lhs dans le residu
%        - Respart_aa:  cellarray (taille Qa x Qa, contenant des matrices 
%            de taille N x N) des termes d'interaction lhs-lhs dans le residu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nombre de termes RHS/ LHS
% --------------------------
Qa = length(AA_decomp);
Ql = length(LL_decomp);
[NbNodes,N] = size(PP);

% calcul des representants de Riesz
% ---------------------------------
BAP = cell(1,Qa);
for i=1:Qa
    BAP{i} = zeros(NbNodes,N);
    for j=1:N
        BAP{i}(:,j) = BB \ (AA_decomp{i}*PP(:,j));
    end
end

% calcul des interactions l-a
% ---------------------------
Respart_la = cell(Ql,Qa);
for q=1:Ql
    for k=1:Qa
        Respart_la{q,k} = BAP{k}'*LL_decomp{q};
    end
end

% calcul des interactions a-a
% ---------------------------
Respart_aa = cell(Qa,Qa);
for q=1:Qa
    for k=q:Qa
        Respart_aa{q,k} = PP'*AA_decomp{k}'*BAP{q};
        % on complete par symetrie:
        Respart_aa{k,q} = Respart_aa{q,k}';
    end
end

end

