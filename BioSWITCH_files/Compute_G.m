function [G,out] = Compute_G(ESP,COM,KNS)


out = Deficiency_computations(ESP,COM,KNS);



lhs = out.A_kPsi; 
rhs = out.rhs; 
Psi = out.Psi; 
Psi_c = out.Psi_c;

alpha_vec = out.alpha_vec;

ESP = out.esp;
COM = out.com;
KNS = out.kns;

n_com = out.n_com;
n_link = out.n_link;
n_esp = out.n_esp;

matrix_N = matrixN(ESP,COM,KNS);

BB = full(localsemiposalg(matrix_N))';

[reduced_lhs, reduced_rhs] = Extract_Lin_Indep_Eqs(lhs,rhs,Psi);

for ii=1:1:n_com-n_link
H(ii)=reduced_lhs(ii)-reduced_rhs(ii);
H_s(ii) = subs(H(ii),Psi,Psi_c);
end

for ii=1:1:n_esp
    eval(sprintf('syms x_%d',ii));
    eval(sprintf('c(%d)=x_%d;',ii,ii));
end

W = BB'*c.';

J1 = jacobian(H_s,c);
J2 = jacobian(H_s,alpha_vec);
J3 = jacobian(W,c);
J4 = jacobian(W,alpha_vec);

G = [J1 J2; J3 J4];
end

% PSY: functions 'localsemiposalg' and 'localgcd' taken from 'sbioconsmoiety'
% -------------------------------------------------------------------------
% Semipositive algorithm. Takes a stoich matrix Stoich.
function T = localsemiposalg(Stoich)
   
   [n, r] = size(Stoich);
   T     = [ speye(n) sparse(Stoich) ]; % zeroth tableau
   
   for j = n+r:-1:n+1 % loop backwards over cols of Stoich
      
      Tnew    = T( T(:,j)==0, 1:j-1 ); % new T: start w/'zero rows' of previous T
      posinds = find( T(:,j) > 0 );    % col vec 
      neginds = find( T(:,j) < 0 );    % col vec
      lni     = length(neginds);
      
      for i = posinds' % loop over positive entries in T(:,j)
         % OR current pos row against all neg rows; vectorize for speed
         zerosets = ~( repmat(T(i,1:n),lni,1) | T(neginds,1:n) ); 
         for k = 1:lni % loop over negative rows
             flags = any( T(:,zerosets(k,:)), 2); % 'incomplete' logical indexing
             flags([i;neginds(k)]) = true;        % ignore these two rows in decision
             if flags
                newrow = -T(neginds(k),j)*T(i,1:j-1) + T(i,j)*T(neginds(k),1:j-1);
                Tnew     = [ Tnew; newrow ]; %#ok<AGROW>
             end
         end
      end % for i = posinds'
      
      T = Tnew;
   end % for j = n+1:n+r
   
   g = zeros(size(T,1),1);
   for c = 1:length(g)
       g(c) = localgcd(T(c,:));
   end
   T = bsxfun(@rdivide,T,g);
   
end % function localsemiposalg

% -------------------------------------------------------------------------
% Find the gcd of a set of nums. Could be cleverer, but this is cheap
% compared to the semipos alg. 
function g = localgcd(v)
   
   if ~isequal(v, round(v))
      g = 1;
      return;
   end

   x = v(v~=0); % consider nonzero elements
   g = 0;
   for i = 1:numel(x)
      g = gcd(g, x(i));
      if g == 1
         break;
      end
   end
end % function localgcd