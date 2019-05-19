function [fout1, fout2] = BioSWITCH_mkfiles_lpsearch_MassCon(network_file)
%[fout1, fout2] = BioSWITCH_mkfiles_lpsearch_MassCon(network_file)
% Write objective function script for deficiency approach analysis
% Write default script to call the optimization solver 


%eval(sprintf('[ESP, COM, KNS] = %s;',network_file));
network_file_function = str2func(network_file);
[ESP, COM, KNS] = network_file_function();

n_esp = size(ESP,2);
matrix_N = matrixN(ESP,COM,KNS);
rankN = rank(matrix_N);
n_cons_laws = n_esp - double(rankN);
n_kns = size(KNS,2);

[G,out] = Compute_G(ESP,COM,KNS);


if out.def<1
    error('Deficiency zero, this approach is only valid for deficiency greater than zero')
end

[k_extend, x_extend, set1, Str_x_extend] = Solve_x_extend(ESP,COM,KNS,network_file);

x_var = x_extend(set1);

if size(x_var,2)>size(x_var,1)
    x_var = x_var.';
end

name_file = sprintf('%s_ObjF_DEF.m',network_file);


fidOut=fopen(name_file,'w');

fprintf(fidOut,'function [f,g]=%s_ObjF_DEF(d_var)\n',network_file);

fprintf(fidOut,'\n');


for ii=1:1:size(k_extend,1)
    fprintf(fidOut,'%s=d_var(%d);\n', char(k_extend(ii)),ii);
end

fprintf(fidOut,'\n');

eval(sprintf('nb_sol=max(size(Str_x_extend.%s));',char(x_var(1))));

if nb_sol> 1
   sprintf('The symbolic system has %d solutions, choose the solution from 1 to %d', nb_sol, nb_sol)
   csol = input('chosen solution:');
else
   csol= 1;
end

for ii=1:1:size(x_var,1)       
    
    eval(sprintf('vec_str{%d}=Str_x_extend.%s(csol);',ii,char(x_var(ii))));
        
    fprintf(fidOut,'%s=%s;\n', char(x_var(ii)),char(vec_str{ii}));
end

fprintf(fidOut,'\n');

fprintf(fidOut,'isn = find(isnan(d_var));\n');
fprintf(fidOut,'isi = find(isinf(d_var));\n');

fprintf(fidOut,'\n');

fprintf(fidOut,'if ~isempty(isn) || ~isempty(isi)\n');
fprintf(fidOut,'\t\t\t f=1e20;\n');
fprintf(fidOut,'\t\t\t g=zeros(%d,1);\n',size(x_var,1)); % check
fprintf(fidOut,'\t\t\t return\n');
fprintf(fidOut,'end\n');

fprintf(fidOut,'\n');

fprintf(fidOut,'G=[\t');
for jj=1:1:size(G,1)
    if jj ~=1
        fprintf(fidOut,'\t');
    end
    for ii=1:1:size(G,1)
        fprintf(fidOut,'%s\t',char(G(jj,ii)));
    end
    if jj ~=size(G,1)
        fprintf(fidOut,'\n');
    end
end
fprintf(fidOut,'];\n');

fprintf(fidOut,'\n');

fprintf(fidOut,'f=1e10*det(G)^2;\n');  % TODO: make the multiplier, here fixed to 1e10, be tunable

fprintf(fidOut,'\n');

cont=0;
for ii=1:1:size(x_var,1)
    if isempty(strfind(char(x_var(ii)),'alpha'))
        cont=cont+1;
    fprintf(fidOut,'g(%d)=%s;\n',cont,char(x_var(ii)));
    end
end

fout1 = fclose(fidOut);



name_file = sprintf('%s_input_lpsearch_MassCon_default.m',network_file);

fidOut2=fopen(name_file,'w');


fprintf(fidOut2,'%% Default Input File for Optimization (networks with mass conservation)\n');

fprintf(fidOut2,'\n');

fprintf(fidOut2,'%% Default Decision Variables\n');
for ii=1:1:size(k_extend,1)    
    fprintf(fidOut2,'%% %s\n',char(k_extend(ii)));    
end

fprintf(fidOut2,'\n');

fprintf(fidOut2,'%% Objective Function File\n');
fprintf(fidOut2,'problem.f=''%s_ObjF_DEF'';\n',network_file);

fprintf(fidOut2,'\n');

fprintf(fidOut2,'%% lower and upper bounds for the decission variables \n');

for ii=1:1:size(k_extend,1)
    k_extend_char{ii}=char(k_extend(ii));
    if k_extend_char{ii}(1)=='k'
        x_L(ii)=1e-2;
        x_U(ii)=1e2;
    elseif k_extend_char{ii}(1)=='x'
        x_L(ii)=1e-2;
        x_U(ii)=1e2;
    else
        x_L(ii)=-1e3;
        x_U(ii)=1e3;
    end
end

for ii=1:1:size(k_extend,1)    
        x_0(ii)=1;         
end

cont=0;
for ii=1:1:size(x_var,1) 
    if isempty(strfind(char(x_var(ii)),'alpha'))
        cont=cont+1;
        c_L(cont)=1e-2;
        c_U(cont)=1e2;  
    end
end      


fprintf(fidOut2,'problem.x_L=[');
for ii=1:1:size(k_extend,1)
     if ii==size(k_extend,1)
        fprintf(fidOut2,'\t %0.2g];\n',x_L(ii));
     else
    fprintf(fidOut2,'\t %0.2g',x_L(ii)); 
     end
end

fprintf(fidOut2,'problem.x_U=[');
for ii=1:1:size(k_extend,1)
     if ii==size(k_extend,1)
        fprintf(fidOut2,'\t %0.2g];\n',x_U(ii));
     else
    fprintf(fidOut2,'\t %0.2g',x_U(ii)); 
     end
end

fprintf(fidOut2,'problem.x_0=[');
for ii=1:1:size(k_extend,1)
     if ii==size(k_extend,1)
        fprintf(fidOut2,'\t %0.2g];\n',x_0(ii));
     else
    fprintf(fidOut2,'\t %0.2g',x_0(ii)); 
     end
end

fprintf(fidOut,'\n');

fprintf(fidOut2,'problem.c_L=[');
for ii=1:1:cont
     if ii==cont
        fprintf(fidOut2,'\t %0.2g];\n',c_L(ii));
     else
    fprintf(fidOut2,'\t %0.2g',c_L(ii)); 
     end
end

fprintf(fidOut2,'problem.c_U=[');
for ii=1:1:cont
     if ii==cont
        fprintf(fidOut2,'\t %0.2g];\n',c_U(ii));
     else
    fprintf(fidOut2,'\t %0.2g',c_U(ii)); 
     end
end

fprintf(fidOut2,'\n');

fprintf(fidOut2,'opts.maxtime=150;\n');
%fprintf(fidOut2,'opts.local.n1=2;\n');
%fprintf(fidOut2,'opts.local.n1=3;\n');

fout2=fclose(fidOut2);


name_file2 = sprintf('%s_param_MassCon.m',network_file);


fidOut3 = fopen(name_file2,'w');

x_var = x_extend(set1);

if size(x_var,2)>size(x_var,1)
    x_var = x_var.';
end

fprintf(fidOut3,'function [x_ss, par_cont] = %s_param_MassCon(d_var)\n',network_file);


fprintf(fidOut3,'\n');

fprintf(fidOut3,'\n');

for ii=1:1:size(k_extend,1)
    fprintf(fidOut3,'%s=d_var(%d);\n', char(k_extend(ii)),ii);
end

fprintf(fidOut3,'\n');

for ii=1:1:size(x_var,1)
    eval(sprintf('vec_str{%d}=Str_x_extend.%s(csol);',ii,char(x_var(ii))));
    fprintf(fidOut3,'%s=%s;\n', char(x_var(ii)),char(vec_str{ii}));
end

fprintf(fidOut3,'\n');

B =  Compute_B(ESP, COM, KNS);

if ~isempty(B)
    
% check version
if verLessThan('matlab','9.2')
    for ii=1:1:n_esp
        Esym(ii) = sym(ESP{ii});
    end
else
    for ii=1:1:n_esp
        Esym(ii) = str2sym(ESP{ii});
    end
end

for ii=1:1:size(B,2)
b(ii) = B(:,ii)'*Esym.';
end

else
    error('Check the network, since there are no mass conservation laws')
end

% check version
if verLessThan('matlab','9.2')
    for ii=1:1:n_esp
        x_ESP(ii)=sym(ESP{ii});
    end
else
    for ii=1:1:n_esp
        x_ESP(ii)=str2sym(ESP{ii});
    end
end

fprintf(fidOut3,'\n');

for ii=1:1:n_esp
    eval(sprintf('syms x_%d',ii)); 
    eval(sprintf('x(%d)=x_%d;',ii,ii));
end

b_x = subs(b,x_ESP,x);

for ii=1:1:n_cons_laws
    fprintf(fidOut3,'C%d=%s;\n',ii,char(b_x(ii)));
end

fprintf(fidOut3,'\n');

for ii=1:1:n_kns
    fprintf(fidOut3,'par_cont(%d)=%s;\n',ii,KNS{ii});
end
for ii=1:1:n_cons_laws
    fprintf(fidOut3,'par_cont(%d)=C%d;\n',ii+n_kns,ii);
end

fprintf(fidOut3,'\n');

for ii=1:1:n_esp
    fprintf(fidOut3,'x_ss(%d)=x_%d;\n',ii,ii);
end

% fprintf(fidOut3,'\n');
% 
% for ii=1:1:n_esp-n_cons_laws
%     fprintf(fidOut3,'idx(%d)=%d;\n',ii,idx(ii));
% end

fclose(fidOut3);
%


%



%
