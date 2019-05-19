function YY = matrixY(ESP,COM)

n_esp=size(ESP,2);
n_com=size(COM,2);

YY=zeros(n_esp,n_com); 


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

end