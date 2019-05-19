function A_k = matrixAk(COM,KNS)

%--------------------------------------------------------------------------
n_com=size(COM,2);
n_rec=size(KNS,2); % new   


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

% COMPUTING MATRIX Ak and KERNEL 
%--------------------------------------------------------------------------
A_k=RR-diag(RR.'*ones(n_com,1));
nucleo=null(A_k);
n_link=size(nucleo,2);  % t=l