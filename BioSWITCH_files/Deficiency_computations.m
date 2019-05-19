function out = Deficiency_computations(ESP,COM,KNS)

% CRNT analysis
% Provided one strong terminal linkage class per linkage class
%--------------------------------------------------------------------------
% Given:
% 1) The vector of species
% 2) The vector of complexes
% 3) The vector of constants
% Compute:
% 1) The deficiency of the network [out.def];
% 2) The number of linkage classes [out.n_link];
% 3) The vector of mass action monomials [out.Psi_c];
% 4) The ODE equations describing the dynamics of the system [out.modeleq];
% 5) The stoichiometric matrix [out.N];
% 6) The molecularity matrix [out.Y];
% 7) The basis of the deficiency subspace [out.basis];
% 8) The matrix A_k
% 9) A_k times Psi (lhs of the linear eq syst) [out.A_kPsi]
%10) alpha*w (rhs of the linear eq syst) [out.rhs]



% COMPUTING MATRIX Y
%--------------------------------------------------------------------------
n_esp=size(ESP,2);
n_com=size(COM,2);
n_rec=size(KNS,2); % new   

YY=zeros(n_esp,n_com); % new


for ii=1:n_esp   
  for jj=1:n_com
           for kk=1:n_esp          
          if kk==ii
              eval(sprintf('%s=1;',ESP{kk}))
          else
              eval(sprintf('%s=0;',ESP{kk}))
          end
           end     
       YY(ii,jj)=eval(COM{jj});
   end
end
r=rank(YY);

% COMPUTING VECTOR Psi
%--------------------------------------------------------------------------
Psi_c=sym(zeros(n_com,1));
esp_v=sym(zeros(n_esp,1));

for ii=1:n_esp
    eval(sprintf('syms x_%d; esp_v(%d)= x_%d;',ii,ii,ii));
    eval(sprintf('x(%d)=x_%d;',ii,ii));
end

for jj=1:n_com
    eval(sprintf('Psi_c(%d)=prod(esp_v.^YY(:,%d));',jj,jj));
end

% SIMBOLIC VECTOR psi 
%--------------------------------------------------------------------------
Psi=sym(zeros(n_com,1));
for ii=1:1:n_com
    eval(sprintf('syms Psi_%d;',ii));
    eval(sprintf('Psi(%d)=Psi_%d;',ii,ii));
end

% COMPUTING MATRIX R (adjacency)
%--------------------------------------------------------------------------
RR=sym(zeros(n_com));

% check version
if verLessThan('matlab','9.2')
    for ii=1:n_rec
        for jj=1:n_com
            for kk=1:n_com
                a=sprintf('k_%d_%d',jj,kk);
                if strcmp(KNS(ii),a)==1
                    RR(kk,jj)=sym(KNS{ii});
                end
            end
        end
    end
else
    for ii=1:n_rec
        for jj=1:n_com
            for kk=1:n_com
                a=sprintf('k_%d_%d',jj,kk);
                if strcmp(KNS(ii),a)==1
                    RR(kk,jj)=str2sym(KNS{ii});
                end
            end
        end
    end
end

% COMPUTING MATRIX M (Kaltenbach notation)
%--------------------------------------------------------------------------
M=sym(zeros(n_com,n_rec));
for ii=1:n_rec
    for jj=1:n_com
        for kk=1:n_com
        a=sprintf('k_%d_%d',jj,kk);
        if strcmp(KNS(ii),a)==1
            M(kk,ii)=1;
            M(jj,ii)=-1;
        end
       end
    end
end 

% COMPUTING MATRIX Ak and KERNEL 
%--------------------------------------------------------------------------
A_k=RR-diag(RR.'*ones(n_com,1));
nucleo=null(A_k);
n_link=size(nucleo,2);  % t=l

% KERNEL REORDERING
%--------------------------------------------------------------------------
% the ordering strategy is the following: first, columns are reordered to
% locate nodes 1, 2 and 3. After that, those columns with ones not located in the
% correct position are multiplied by the corresponding factor.

% reordering columns (ver esto)
%-------------------------------------------------------------------------
nucleo_reord=nucleo;
for ii=1:1:n_link
    for jj=1:1:n_link
        if [nucleo(ii,jj) == 0] == 0
            nucleo_reord(:,ii)=nucleo(:,jj);
        end
    end
end

% factorizing
%-------------------------------------------------------------------------
for ii=1:1:n_link 
    if ([nucleo_reord(ii,ii)-1 == 0] == 0 && nucleo_reord(ii,ii) ~= 0 )
        nucleo_reord(:,ii)=nucleo_reord(:,ii)/nucleo_reord(ii,ii);
    end
end

% BASIS for D_delta
% -------------------------------------------------------------------------
M=double(M);
basis=null([YY; null(M.','r').'],'r');


% DIMENSION OF THE STOICHIOMETRIC SUBSPACE
%--------------------------------------------------------------------------
rank_YAk = double(rank(YY*A_k));

% DEFICIENCY
%--------------------------------------------------------------------------
def = size(A_k,1)-size(nucleo,2)-rank_YAk;

% 

if def>0
    for ii=1:1:def
        eval(sprintf('syms alpha_%d \n', ii));
    end
    
    
    for ii=1:1:def
        eval(sprintf('alpha_vec(%d)=alpha_%d; \n', ii,ii));
    end
    
           
    if size(alpha_vec,1)<size(alpha_vec,2)
        alpha_vec = alpha_vec.';
    end

    if max(size(basis(1,:))) ~= max(size(alpha_vec))
        disp('The network is not uniterminal, we recommend the injectivity based analysis')   
        out.def = def;
        out.n_link =n_link;
        out.n_com=n_com;
        out.n_esp=n_esp;
        out.n_rec=n_rec;
        out.Psi_c = Psi_c;
        out.Psi = Psi;
        out.modeleq = YY*A_k*Psi_c;
        out.N = YY*M;
        out.esp=ESP;
        out.com=COM;
        out.kns=KNS;
        out.Y=YY;      
        out.A_k=A_k;
        out.A_kPsi = A_k*Psi;      
        
        return
    end
    
    for ii=1:1:n_com       
        eval(sprintf('rhs(%d) = basis(%d,:)*alpha_vec;',ii,ii));
    end
    
     if size(rhs,1)<size(rhs,2)
        rhs = rhs.';
    end

else
    for ii=1:1:n_com
        eval(sprintf('rhs = zeros(%d,1);',n_com));
    end
    alpha_vec=[];
end


out.def = def;
out.n_link =n_link;
out.n_com=n_com;
out.n_esp=n_esp;
out.n_rec=n_rec;
out.Psi_c = Psi_c;
out.Psi = Psi;
out.modeleq = YY*A_k*Psi_c;
out.N = YY*M;
out.esp=ESP;
out.com=COM;
out.kns=KNS;
out.Y=YY;
out.basis=basis;
out.A_k=A_k;
out.A_kPsi = A_k*Psi;
out.rhs = rhs;
out.alpha_vec = alpha_vec;




