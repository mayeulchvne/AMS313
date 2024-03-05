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

% calcul des representants de Riesz
% ---------------------------------


% calcul des interactions l-a
% ---------------------------
Respart_la = cell(Ql,Qa);


% calcul des interactions a-a
% ---------------------------
Respart_aa = cell(Qa,Qa);

error('RB_compute_respart_la_aa() not yet implemented');

end

