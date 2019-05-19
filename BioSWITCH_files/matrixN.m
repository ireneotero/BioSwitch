function matrixN = matrixN(ESP,COM,KNS)
 % function matrixN = matrixN(ESP,COM,KNS)
 % computes stoichoimetric matrix N 

% COMPUTING MATRIX Y
%--------------------------------------------------------------------------
n_esp=size(ESP,2);
n_com=size(COM,2);
n_rec=size(KNS,2);   


YY=zeros(n_esp,n_com); 
for ii=1:n_esp   
  for jj=1:n_com      
           for kk=1:n_esp          
           if kk==ii
              eval(sprintf('%s=1;',ESP{kk}));
           else
              eval(sprintf('%s=0;',ESP{kk}));
           end
           end    
       YY(ii,jj)=eval(COM{jj});
   end
end


% COMPUTING MATRIX M (Kaltenbach notation)  
%--------------------------------------------------------------------------
M=zeros(n_com,n_rec);
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


matrixN = YY*M;

