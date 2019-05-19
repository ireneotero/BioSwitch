function [fout1, fout2] = BioSWITCH_mkfiles_lpsearch_SemiDiff(network_file,input_species_vector)
%function fout = BioSWITCH_mkfiles_lpsearch_SemiDiff(network_file,input_species_vector)
%
% Write objective function script for injectivity approach analysis
% Write default script to call the optimization solver 
%
% network_file: name of the script containing the closed network
%  input_species_vector: row vector containing indices of input species (one per c.law). The indices must be in the right
%  order according to the labeling of the species. For example:
% if the species are ordered such that { 'A', 'B', 'C', 'D' }
% and the conservation laws are
%  C1 = A + B + C
%  C2 = 2A + D 
% If we have in the input C and A, the species vector is [1,3]

%eval(sprintf('[ESP, COM, KNS] = %s;',network_file));
network_file_function = str2func(network_file);
[ESP, COM, KNS] = network_file_function();

n_esp=size(ESP,2);


name_file = sprintf('%s_ObjF_INJ.m',network_file);

fidOut = fopen(name_file,'w');

fprintf(fidOut,'function [f,g]=%s_ObjF_INJ(d_var)\n',network_file);

fprintf(fidOut,'\n');

fprintf(fidOut,'mu=d_var;\n');

fprintf(fidOut,'\n');

Y = matrixY(ESP,COM);

fprintf(fidOut,'Y=[');
for ii = 1:1:size(Y,1)-1
    for jj=1:1:size(Y,2)
        if jj==size(Y,2)
            fprintf(fidOut,'\t %d\n',Y(ii,jj));
        else
            fprintf(fidOut,'\t %d',Y(ii,jj));
        end
    end
end
for jj = 1:1:size(Y,2)
    if jj==size(Y,2)
        fprintf(fidOut,'\t %d];\n',Y(size(Y,1),jj));
    else
        fprintf(fidOut,'\t %d',Y(size(Y,1),jj));
    end
end

[Yr, Nr] = matrixYrNr(ESP,COM,KNS,input_species_vector);


fprintf(fidOut,'\n');
fprintf(fidOut,'Yr=[');
for ii = 1:1:size(Yr,1)-1
    for jj = 1:1:size(Yr,2)
        if jj==size(Yr,2)
            fprintf(fidOut,'\t %d\n',Yr(ii,jj));
        else
            fprintf(fidOut,'\t %d',Yr(ii,jj));
        end
    end
end
for jj = 1:1:size(Yr,2)
    if jj==size(Yr,2)
        fprintf(fidOut,'\t %d];\n',Yr(size(Yr,1),jj));
    else
        fprintf(fidOut,'\t %d',Yr(size(Yr,1),jj));
    end
end

fprintf(fidOut,'\n');

fprintf(fidOut,'Nr=[');
for ii = 1:1:size(Nr,1)-1
    for jj = 1:1:size(Nr,2)
        if jj==size(Nr,2)
            fprintf(fidOut,'\t %d\n',Nr(ii,jj));
        else
            fprintf(fidOut,'\t %d',Nr(ii,jj));
        end
    end
end
for jj = 1:1:size(Nr,2)
    if jj==size(Nr,2)
        fprintf(fidOut,'\t %d];\n',Nr(size(Nr,1),jj));
    else
        fprintf(fidOut,'\t %d',Nr(size(Nr,1),jj));
    end
end

fprintf(fidOut,'\n');

fprintf(fidOut,'J = Nr*diag(mu)*Yr'';\n');

fprintf(fidOut,'\n');

fprintf(fidOut,'f=1e12*det(J)^2;\n');  %TODO 1e12 tunable

fprintf(fidOut,'\n');
 
fprintf(fidOut,'p = -Nr*mu'';\n');
 
fprintf(fidOut,'\n');
 
fprintf(fidOut,'%% equality constraints\n');

no_input_species_vector=setxor(1:n_esp,input_species_vector); 

if ~isempty(no_input_species_vector)
for ii=1:1:max(size(no_input_species_vector))
    fprintf(fidOut,'g(%d)=p(%d);\n',ii,no_input_species_vector(ii));
end
end

fprintf(fidOut,'\n');

fprintf(fidOut,'%% inequality constraints\n');
for ii=1:1:max(size(input_species_vector))
    fprintf(fidOut,'g(%d)=p(%d);\n',ii+max(size(no_input_species_vector)),input_species_vector(ii));
end


fout1 = fclose(fidOut);


name_file = sprintf('%s_input_lpsearch_SemiDiff_default.m',network_file);

fidOut2=fopen(name_file,'w');


fprintf(fidOut2,'%% Default Input File for Optimization (semidiffusive case)\n');

fprintf(fidOut2,'\n');
 
n_dvar=size(Yr,2);
for ii=1:1:max(size(input_species_vector))
    fprintf(fidOut2,'%% input_species %d = %s (index %d) \n', ii, ESP{input_species_vector(ii)}, input_species_vector(ii));
end

fprintf(fidOut2,'\n');
fprintf(fidOut2,'%% Default Decision Variables\n');
for ii=1:1:n_dvar    
     fprintf(fidOut2,'%% mu_%d\n',ii);    
end


fprintf(fidOut2,'\n');

fprintf(fidOut2,'%% Objective Function File\n');
fprintf(fidOut2,'problem.f=''%s_ObjF_INJ'';\n',network_file);
 
fprintf(fidOut2,'\n');
 
fprintf(fidOut2,'%% lower and upper bounds for the decision variables \n');

for ii=1:1:n_dvar
         x_L(ii)=1e-2;
         x_U(ii)=1e2;
end

for ii=1:1:n_dvar
         x_0(ii)=1;
end

fprintf(fidOut2,'problem.x_L=[');
for ii=1:1:n_dvar
      if ii==n_dvar
         fprintf(fidOut2,'\t %0.2g];\n',x_L(ii));
      else
     fprintf(fidOut2,'\t %0.2g',x_L(ii)); 
      end
 end

fprintf(fidOut2,'problem.x_U=[');
for ii=1:1:n_dvar
      if ii==n_dvar
         fprintf(fidOut2,'\t %0.2g];\n',x_U(ii));
      else
     fprintf(fidOut2,'\t %0.2g',x_U(ii)); 
      end     
end


fprintf(fidOut2,'problem.x_0=[');
 for ii=1:1:n_dvar
      if ii==n_dvar
         fprintf(fidOut2,'\t %0.2g];\n',x_0(ii));
      else
     fprintf(fidOut2,'\t %0.2g',x_0(ii)); 
      end
 end
 
 
fprintf(fidOut,'\n');

fprintf(fidOut2,'%% lower and upper bounds for inequality constraints \n');

for ii=1:1:max(size(input_species_vector))
         c_L(ii)=0;
         c_U(ii)=Inf;
end

fprintf(fidOut2,'problem.c_L=[');
for ii=1:1:max(size(input_species_vector))
     if ii==max(size(input_species_vector))
        fprintf(fidOut2,'\t %0.2g];\n',c_L(ii));
     else
    fprintf(fidOut2,'\t %0.2g',c_L(ii)); 
     end
end

fprintf(fidOut2,'problem.c_U=[');
for ii=1:1:max(size(input_species_vector))
     if ii==max(size(input_species_vector))
        fprintf(fidOut2,'\t %0.2g];\n',c_U(ii));
     else
    fprintf(fidOut2,'\t %0.2g',c_U(ii)); 
     end
end


fprintf(fidOut,'\n');

fprintf(fidOut2,'problem.neq=%d;\n',n_esp-max(size((input_species_vector))));

fprintf(fidOut,'\n');

fprintf(fidOut2,'opts.maxtime=150;\n');
%fprintf(fidOut2,'opts.local.n1=2;\n');
%fprintf(fidOut2,'opts.local.n1=3;\n');

fprintf(fidOut,'\n');

fout2=fclose(fidOut2);



