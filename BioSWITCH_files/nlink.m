function n_link = nlink(ESP,COM,KNS)
% function n_link = nlink(ESP,COM,KNS)
% computes the number of linkage classes
%

n_esp=size(ESP,2);
n_com=size(COM,2);
n_rec=size(KNS,2); 


for ii=1:n_esp   
  for jj=1:n_com
           for kk=1:n_esp          
          if kk==ii
              eval(sprintf('%s=1;',ESP{kk}))
          else
              eval(sprintf('%s=0;',ESP{kk}))
          end
           end      
   end
end


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
                    RR(kk,jj)=sym(KNS{ii});
                end
            end
        end
    end
end
    

A_k=RR-diag(RR.'*ones(n_com,1));
nucleo=null(A_k);
n_link=size(nucleo,2);  % t=l

