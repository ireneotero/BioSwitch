function B = Compute_B(ESP,COM,KNS)

matrix_N = matrixN(ESP,COM,KNS);

B = full(localsemiposalg(matrix_N))';
end

% functions 'localsemiposalg' and 'localgcd' are taken from 'sbioconsmoiety', which is part of SimBiology
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