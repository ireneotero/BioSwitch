function [fout1, fout2] = BioSWITCH_mkfiles_ncontin_SemiDiff(network_file,input_species_vector)
%  [fout1, fout2] = BioSWITCH_mkfiles_ncontin_SemiDiff(network_file,input_species_vector)
%  network_file: name of the script containing the network (only true
%  reactions)

%  input_species_vector: vector with indices of species in the input.
%  The indices must be in the right
%  order according to the labeling of the species from low to high
%  indicies.
%  For example:
%  if the species are ordered such that { 'A', 'B', 'C', 'D' }
%  and the conservation laws are
%  C1 = A + B + C
%  C2 = 2A + D
%  If we want to introduce in the input C and A, the input_species_vector is [1,3]

%eval(sprintf('[ESP, COM, KNS] = %s;',network_file));
network_file_function = str2func(network_file);
[ESP, COM, KNS] = network_file_function();

%% first file (ODE file for continuation with Cl-Matcont)
name_file = sprintf('%s_ode_SemiDiff.m',network_file);

fidOut3 = fopen(name_file,'w');

n_esp = size(ESP,2);
n_kns = size(KNS,2);


n_esp_in = max(size(input_species_vector));


fprintf(fidOut3,'function out=%s_ode_SemiDiff\n',network_file);

fprintf(fidOut3,'\n');

fprintf(fidOut3,'out{1}=[];\n');
fprintf(fidOut3,'out{2}=@fun_eval;\n');
fprintf(fidOut3,'out{3}=@jacobian;\n');
fprintf(fidOut3,'out{4}=[];\n');
fprintf(fidOut3,'out{5}=[];\n');
fprintf(fidOut3,'out{6}=[];\n');
fprintf(fidOut3,'out{7}=[];\n');
fprintf(fidOut3,'out{8}=[];\n');
fprintf(fidOut3,'out{9}=[];\n');

fprintf(fidOut3,'\n');

fprintf(fidOut3,'function dxdt = fun_eval(t,x,');

for ii=1:1:n_kns+n_esp+n_esp_in-1
    fprintf(fidOut3,'par%d,',ii);
end
fprintf(fidOut3,'par%d)\n\n',n_kns+n_esp+n_esp_in);


out = Deficiency_computations(ESP, COM, KNS);

modeleq = out.modeleq;


for ii=1:1:n_esp
    fprintf(fidOut3,'%s=x(%d);\n',ESP{ii},ii);
end

fprintf(fidOut3,'\n');

for ii=1:1:n_kns
    fprintf(fidOut3,'%s=par%d;\n',KNS{ii},ii);
end

for ii=1:1:n_esp
    fprintf(fidOut3,'kdeg_%s=par%d;\n',ESP{ii},n_kns+ii);
end

for ii=1:1:n_esp_in
    fprintf(fidOut3,'K%d=par%d;\n',input_species_vector(ii),n_kns+n_esp+ii);
end

fprintf(fidOut3,'\n');



fprintf(fidOut3,'\n');

for ii=1:1:n_esp
    eval(sprintf('syms x_%d',ii));
    eval(sprintf('x(%d)=x_%d;',ii,ii));
end

% check version
if verLessThan('matlab','9.2')
    modeleq2 = subs(modeleq,x,sym(ESP));
else
    modeleq2 = subs(modeleq,x,str2sym(ESP));
end

for ii=1:1:n_esp
    eval(sprintf('syms %s',ESP{ii}));
    eval(sprintf('syms kdeg_%s',ESP{ii}));
    
    eval(sprintf('modeleq3(ii)= modeleq2(ii) - kdeg_%s*sym(%s);',ESP{ii},ESP{ii}))
end



fprintf(fidOut3,'\n');

for ii=1:1:n_esp
    if any(input_species_vector==ii)
        fprintf(fidOut3,'d%sdt = K%d + %s -kdeg_%s*%s;\n',ESP{ii},ii,char(modeleq2(ii)),ESP{ii},ESP{ii});
    else
        fprintf(fidOut3,'d%sdt = %s -kdeg_%s*%s;\n',ESP{ii},char(modeleq2(ii)),ESP{ii},ESP{ii});
    end
end

fprintf(fidOut3,'\n');

fprintf(fidOut3,'dxdt=[')
for ii=1:1:n_esp-1
    fprintf(fidOut3,'d%sdt; ',ESP{ii});
end
fprintf(fidOut3,'d%sdt];\n',ESP{n_esp});clc

for ii=1:1:n_esp
    eval(sprintf('syms x_%d',ii));
    eval(sprintf('c(%d)=x_%d;',ii,ii));
end


fprintf(fidOut3,'\n');

fprintf(fidOut3,'function [dfdx, detJAC, eigJAC] = jacobian(t,x,');

for ii=1:1:n_kns+n_esp+n_esp_in-1
    fprintf(fidOut3,'par%d,',ii);
end
fprintf(fidOut3,'par%d)\n\n',n_kns+n_esp+n_esp_in);


for ii=1:1:n_esp
    x_ESP(ii)=sym(ESP{ii});
end

fprintf(fidOut3,'\n');

for ii=1:1:n_esp
    fprintf(fidOut3,'%s=x(%d);\n',ESP{ii},ii);
end

fprintf(fidOut3,'\n');
JAC = jacobian(modeleq3,x_ESP);


for ii=1:1:n_kns
    fprintf(fidOut3,'%s=par%d;\n',KNS{ii},ii);
end

for ii=1:1:n_esp
    fprintf(fidOut3,'kdeg_%s=par%d;\n',ESP{ii},n_kns+ii);
end

for ii=1:1:n_esp_in
     fprintf(fidOut3,'K%d=par%d;\n',input_species_vector(ii),n_kns+n_esp+ii);
end


fprintf(fidOut3,'\n');


fprintf(fidOut3,'\n');

for ii=1:1:n_esp
    eval(sprintf('syms x_%d',ii));
    eval(sprintf('x(%d)=x_%d;',ii,ii));
end

fprintf(fidOut3,'\n');

fprintf(fidOut3,'dfdx=[\t');
for ii=1:1:n_esp
    for jj=1:1:n_esp
        if jj== n_esp && ii~=n_esp
            fprintf(fidOut3,'%s\n', char(JAC(ii,jj)));
        elseif jj==n_esp && ii==n_esp
            fprintf(fidOut3,'%s];\n', char(JAC(ii,jj)));
        else
            fprintf(fidOut3,'%s,\t', char(JAC(ii,jj)));
        end
    end
end
%
%
fprintf(fidOut3,'\n');
%
fprintf(fidOut3,'detJAC = det(dfdx);\n');
%
fprintf(fidOut3,'\n');
%
fprintf(fidOut3,'eigJAC = eig(dfdx);\n');
%

fout1 = fclose(fidOut3);


%% second file (parameters)

name_file2 = sprintf('%s_param_SemiDiff.m',network_file);


fidOut4 = fopen(name_file2,'w');

[Yr, Nr] = matrixYrNr(ESP,COM,KNS,input_species_vector);

fprintf(fidOut4,'function [xss, par_cont, idx] = %s_param_SemiDiff(d_var)\n',network_file);

n_dvar=size(Yr,2);

ESP_v=sym(ESP).';

for ii=1:n_dvar-n_esp
    mu(ii)=KNS{ii}*prod(ESP_v.^Yr(:,ii));
end

for ii=1:n_esp
    eval(sprintf('mu(ii+n_dvar-n_esp)=kdeg_%s*%s;',ESP{ii},ESP{ii}))
end

fprintf(fidOut4,'\n');

for ii=1:1:n_dvar
    fprintf(fidOut4,'mu(%d)=d_var(%d);\n',ii,ii);
end

fprintf(fidOut4,'\n');

for ii=1:1:n_dvar
    fprintf(fidOut4,'%% mu(%d)=%s;\n',ii,char(mu(ii)));
end

fprintf(fidOut4,'\n');

for ii=1:1:n_esp
    fprintf(fidOut4,'xss(%d)=%d;\n',ii,1);
end

fprintf(fidOut4,'\n');

p = -Nr*mu.';

for ii=1:1:n_esp_in
    eval(sprintf('K%d=p(%d);',input_species_vector(ii), input_species_vector(ii)));
end

fprintf(fidOut4,'\n');

for ii=1:1:n_dvar-n_esp
    fprintf(fidOut4,'par_cont(%d)=mu(%d); %% %s\n',ii,ii,KNS{ii});
end
for ii=1:1:n_esp
    fprintf(fidOut4,'par_cont(%d)=mu(%d); %% kdeg_%s\n',ii+n_dvar-n_esp,ii+n_dvar-n_esp,ESP{ii});
end

fprintf(fidOut4,'\n');

for ii=1:1:n_dvar-n_esp
    fprintf(fidOut4,'%s=par_cont(%d); \n',KNS{ii},ii);
end
for ii=1:1:n_esp
    fprintf(fidOut4,'kdeg_%s=par_cont(%d); \n',ESP{ii},ii+n_dvar-n_esp);
end

fprintf(fidOut4,'\n');

for ii=1:1:n_esp
    fprintf(fidOut4,'%s = xss(%d);\n', ESP{ii},ii);
end

fprintf(fidOut4,'\n');
for ii=1:1:n_esp_in
    fprintf(fidOut4,'K%d=%s;\n',input_species_vector(ii), char(p(input_species_vector(ii))));
end

fprintf(fidOut4,'\n');
for ii=1:1:n_esp_in
    fprintf(fidOut4,'par_cont(%d)=K%d;\n',ii+n_dvar,input_species_vector(ii));
end

fprintf(fidOut4,'\n');
for ii=1:1:n_esp
    fprintf(fidOut4,'idx(%d)=%d;\n',ii,ii);
end

fprintf(fidOut4,'\n');


fclose(fidOut4);

%% third file (main file for continuation with Cl-Matcont)

name_file = sprintf('%s_input_cont_SemiDiff_default.m',network_file);
 
fidOut5 = fopen(name_file,'w');
 
n_esp = size(ESP,2);
n_kns = size(KNS,2);
n_esp_in = max(size(input_species_vector));
 
 
fprintf(fidOut5,'%% Default file for Equibrium continuation with CL_Matcont (Semi-Diffusive network)\n');
fprintf(fidOut5,'clear all\n');
fprintf(fidOut5,'clc \n\n');
 
fprintf(fidOut5,'load %s_input_lpsearch_SemiDiff_default_Results.mat\n',network_file);
fprintf(fidOut5,'xbest = Results.xbest;\n');
fprintf(fidOut5,'odefunc = %s_ode_SemiDiff;\n',network_file);
fprintf(fidOut5,'[xss, par_cont,idx] = %s_param_SemiDiff(xbest);\n',network_file);
fprintf(fidOut5,'xss=xss(idx)'';\n');
fprintf(fidOut5,'x0 = xss;\n\n');

fprintf(fidOut5,'%% parameter values\n');
 
for ii=1:1:n_kns+n_esp+n_esp_in
     fprintf(fidOut5,'par%d = par_cont(%d);\n',ii,ii);
end
 
fprintf(fidOut5,'\n');
fprintf(fidOut5,'%% integration\n');
fprintf(fidOut5,'tspan = 0:0.1:1000;\n');
fprintf(fidOut5,'[t,x] = ode45(odefunc{2},tspan, x0,[],');
 
for ii=1:1:n_kns+n_esp+n_esp_in-1
     fprintf(fidOut5,'par%d,',ii);
end
 
fprintf(fidOut5,'par%d);\n\n',n_kns+n_esp+n_esp_in);
fprintf(fidOut5,'%% check steady state condition\n');
fprintf(fidOut5,'figure(1)\n');
fprintf(fidOut5,'plot(t,x)\n\n');

fprintf(fidOut5,'%% continuation with CL_matcont\n');
% 
fprintf(fidOut5,'p=[');
for ii=1:1:n_kns+n_esp+n_esp_in-1
    fprintf(fidOut5,'par%d,',ii);
end
 
fprintf(fidOut5,'par%d];\n\n',n_kns+n_esp+n_esp_in);
fprintf(fidOut5,'%% index of the parameter taken as continuation parameter\n');
fprintf(fidOut5,'ap1=1;\n\n');

fprintf(fidOut5,'%% Initialize\n');
fprintf(fidOut5,'[x0,v0]=init_EP_EP(@%s_ode_SemiDiff,xss,p, ap1);\n\n',network_file);

fprintf(fidOut5,'%% Options of the continuer\n');
fprintf(fidOut5,'opt=contset;\n');
fprintf(fidOut5,'opt=contset(opt,''VarTolerance'',1e-8);\n');
fprintf(fidOut5,'opt=contset(opt,''FunTolerance'',1e-8);\n');
fprintf(fidOut5,'opt=contset(opt,''MaxNumPoints'',300);\n');
fprintf(fidOut5,'opt=contset(opt,''Singularities'',1);\n\n');

fprintf(fidOut5,'%% Forward continuation from xss with respect to ap1\n');
fprintf(fidOut5,'[x,v,s,h,f]=cont(@equilibrium,x0,[],opt);\n\n');

fprintf(fidOut5,'opt=contset(opt,''Backward'',1);\n');
fprintf(fidOut5,'%% Backward continuation from xss with respect to ap1\n');
fprintf(fidOut5,'[xb,vb,sb,hb,fb]=cont(@equilibrium,x0,[],opt);\n');
fprintf(fidOut5,'figure(2)\n\n');
fprintf(fidOut5,'%% plot forward branch\n'); 
fprintf(fidOut5,'plot(x(end,:),x(1,:),''b-'')\n');
fprintf(fidOut5,'hold on\n\n');
fprintf(fidOut5,'%% plot backward branch\n');
fprintf(fidOut5,'plot(xb(end,:),xb(1,:),''r-'')\n\n');
 
fprintf(fidOut5,'%% plot starting point\n');
fprintf(fidOut5,'plot(xb(end,1),xb(1,1),''k*'')\n');
fprintf(fidOut5,'plot(x(end,1),x(1,1),''r*'')\n');

fclose(fidOut5);

