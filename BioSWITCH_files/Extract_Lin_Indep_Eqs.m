% Extract the linearly independent linear equations (in psi)
function [reduced_lhs, reduced_rhs] = Extract_Lin_Indep_Eqs(lhs,rhs,vars)

 warning('off','symbolic:sym:isAlways:TruthUnknown')      
 
% convert the set of linear equations in the symbolic variables 'vars' to
% matrix form (Note: b will always be zero)
[matr_form_lhs,b] = equationsToMatrix(lhs, vars);

% extend the lhs matrix (matr_form_A_kPsi) by the (minus)rhs before employing
% Gaussian elimination
M_extended_lhs_minus_rhs = horzcat(matr_form_lhs, -rhs);

% Reduced row echelon form of matrix (Gauss-Jordan elimination)
reduced_M = rref(M_extended_lhs_minus_rhs);

% keep only non-zero rows
reduced_M( all( isAlways(reduced_M==0) ,2) ,:) = [];

reduced_lhs = reduced_M(:,1:end-1)*vars;

reduced_rhs= reduced_M(:,end);

