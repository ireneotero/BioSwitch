function [fout1] = BioSWITCH_mkfiles_ncontin_MassCon(network_file,species_vector)
%  [fout1, fout2] = BioSWITCH_mkfiles_ncontin_MassCon(network_file,species_vector)

%  network_file: name of the script containing the network, for example 'EX1' 

%  species_vector: vector with indices of species to suppress from the
%  ODEs (one per mass conservation law). The indices must be in the right
%  order according to the labeling of the conservation laws. For example:
% if the species are ordered such that { 'A', 'B', 'C', 'D' }
% and the conservation laws are
%  C1 = A + B + C
%  C2 = 2A + D 
% If we want to supress C and A, the species vector is [3,1]

%eval(sprintf('[ESP, COM, KNS] = %s;',network_file));
network_file_function = str2func(network_file);
[ESP, COM, KNS] = network_file_function();


name_file = sprintf('%s_ode_MassCon.m',network_file);

fidOut3 = fopen(name_file,'w');

n_esp = size(ESP,2);
n_kns = size(KNS,2);
n_cons_laws = max(size(species_vector));

fprintf(fidOut3,'function out=%s_ode_MassCon\n',network_file);

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

for ii=1:1:n_kns+n_cons_laws-1
    fprintf(fidOut3,'par%d,',ii);
end
fprintf(fidOut3,'par%d)\n\n',n_kns+n_cons_laws);

out = Deficiency_computations(ESP, COM, KNS);

modeleq = out.modeleq;

idx = setxor(1:n_esp,species_vector);

for ii=1:1:n_esp-n_cons_laws
    fprintf(fidOut3,'%s=x(%d);\n',ESP{idx(ii)},ii);
end

fprintf(fidOut3,'\n');

for ii=1:1:n_kns
    fprintf(fidOut3,'%s=par%d;\n',KNS{ii},ii);
end

fprintf(fidOut3,'\n');

for ii=1:1:n_cons_laws
    fprintf(fidOut3,'C%d=par%d;\n',ii,ii+n_kns);
end


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


modeleq_red = modeleq2(idx);

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

for ii=1:1:size(B,2)
eq{ii} = sprintf('%s=C%d;\n',char(b(ii)),ii);
end

% check version
if verLessThan('matlab','9.2') 
    for ii=1:1:n_cons_laws
        xx(ii) = solve(sym(eq{ii}),sym(ESP{species_vector(ii)}));
    end
else
    for ii=1:1:n_cons_laws
        xx(ii) = solve(str2sym(eq{ii}),str2sym(ESP{species_vector(ii)}));
    end
end

for ii=1:1:n_cons_laws
fprintf(fidOut3,'%s = %s;\n',ESP{species_vector(ii)},char(xx(ii)));
end

fprintf(fidOut3,'\n');

for ii=1:1:n_esp-n_cons_laws
fprintf(fidOut3,'d%sdt = %s;\n',ESP{idx(ii)},char(modeleq2(idx(ii))));
end

fprintf(fidOut3,'\n');

fprintf(fidOut3,'dxdt=[')
for ii=1:1:n_esp-n_cons_laws-1
    fprintf(fidOut3,'d%sdt; ',ESP{idx(ii)});
end
fprintf(fidOut3,'d%sdt];\n',ESP{idx(n_esp-n_cons_laws)});clc

for ii=1:1:n_esp
    eval(sprintf('syms x_%d',ii)); 
    eval(sprintf('c(%d)=x_%d;',ii,ii));
end


fprintf(fidOut3,'\n');

fprintf(fidOut3,'function [dfdx, detJAC, eigJAC] = jacobian(t,x,');

for ii=1:1:n_kns+n_cons_laws-1
    fprintf(fidOut3,'par%d,',ii);
end
fprintf(fidOut3,'par%d)\n\n',n_kns+n_cons_laws);


for ii=1:1:n_cons_laws
    eval(sprintf('syms C%d',ii));
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

for ii=1:1:n_esp-n_cons_laws
    fprintf(fidOut3,'%s=x(%d);\n',ESP{idx(ii)},ii);
end

fprintf(fidOut3,'\n');
modeleq_red_subs = subs(modeleq_red, x_ESP(species_vector), xx);
JAC = jacobian(modeleq_red_subs,x_ESP(idx));

for ii=1:1:n_kns
    fprintf(fidOut3,'%s=par%d;\n',KNS{ii},ii);
end

fprintf(fidOut3,'\n');

for ii=1:1:n_cons_laws
    fprintf(fidOut3,'C%d=par%d;\n',ii,ii+n_kns);
end


fprintf(fidOut3,'\n');

for ii=1:1:n_esp
    eval(sprintf('syms x_%d',ii)); 
    eval(sprintf('x(%d)=x_%d;',ii,ii));
end

fprintf(fidOut3,'\n');

fprintf(fidOut3,'dfdx=[\t');
 for ii=1:1:n_esp-n_cons_laws
     for jj=1:1:n_esp-n_cons_laws
         if jj== n_esp-n_cons_laws && ii~=n_esp-n_cons_laws
             fprintf(fidOut3,'%s\n', char(JAC(ii,jj)));
         elseif jj==n_esp-n_cons_laws && ii==n_esp-n_cons_laws
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



name_file = sprintf('%s_input_cont_MassCon_default.m',network_file);
fidOut4=fopen(name_file,'w');
fprintf(fidOut4,'%% Default file for Equibrium continuation with CL_Matcont (network with mass conservation)\n');
fprintf(fidOut4,'clear all\n');
fprintf(fidOut4,'clc \n\n');
 
fprintf(fidOut4,'load %s_input_lpsearch_MassCon_default_Results.mat\n',network_file);
fprintf(fidOut4,'xbest = Results.xbest;\n');
fprintf(fidOut4,'odefunc = %s_ode_MassCon;\n',network_file);
fprintf(fidOut4,'[xss, par_cont] = %s_param_MassCon(xbest);\n',network_file);
fprintf(fidOut3,'\n');
 
for ii=1:1:n_esp-n_cons_laws
     fprintf(fidOut4,'idx(%d)=%d;\n',ii,idx(ii));
end

fprintf(fidOut4,'xss=xss(idx)'';\n');
fprintf(fidOut4,'x0 = xss;\n\n');

fprintf(fidOut4,'%% parameter values\n');
 
for ii=1:1:n_kns+n_cons_laws
     fprintf(fidOut4,'par%d = par_cont(%d);\n',ii,ii);
end


 
fprintf(fidOut4,'\n');
fprintf(fidOut4,'%% integration\n');
fprintf(fidOut4,'tspan = 0:0.1:1000;\n');
fprintf(fidOut4,'[t,x] = ode45(odefunc{2},tspan, x0,[],');
 
for ii=1:1:n_kns+n_cons_laws-1
     fprintf(fidOut4,'par%d,',ii);
end
 

fprintf(fidOut4,'par%d);\n\n',n_kns+n_cons_laws);
fprintf(fidOut4,'%% check steady state condition\n');
fprintf(fidOut4,'figure(1)\n');
fprintf(fidOut4,'plot(t,x)\n\n');

fprintf(fidOut4,'%% continuation with CL_matcont\n');
% 
fprintf(fidOut4,'p=[');
for ii=1:1:n_kns+n_cons_laws-1
    fprintf(fidOut4,'par%d,',ii);
end
fprintf(fidOut4,'par%d];\n\n',n_kns+n_cons_laws);

fprintf(fidOut4,'%% index of the parameter taken as continuation parameter\n');
fprintf(fidOut4,'ap1=%d;\n\n',n_kns+n_cons_laws);

fprintf(fidOut4,'%% Initialize\n');
fprintf(fidOut4,'[x0,v0]=init_EP_EP(@%s_ode_MassCon,xss,p, ap1);\n\n',network_file);

fprintf(fidOut4,'%% Options of the continuer\n');
fprintf(fidOut4,'opt=contset;\n');
fprintf(fidOut4,'opt=contset(opt,''VarTolerance'',1e-8);\n');
fprintf(fidOut4,'opt=contset(opt,''FunTolerance'',1e-8);\n');
fprintf(fidOut4,'opt=contset(opt,''MaxNumPoints'',300);\n');
fprintf(fidOut4,'opt=contset(opt,''Singularities'',1);\n\n');

fprintf(fidOut4,'%% Forward continuation from xss with respect to ap1\n');
fprintf(fidOut4,'[x,v,s,h,f]=cont(@equilibrium,x0,[],opt);\n\n');

fprintf(fidOut4,'opt=contset(opt,''Backward'',1);\n');
fprintf(fidOut4,'%% Backward continuation from xss with respect to ap1\n');
fprintf(fidOut4,'[xb,vb,sb,hb,fb]=cont(@equilibrium,x0,[],opt);\n');
fprintf(fidOut4,'figure(2)\n\n');
fprintf(fidOut4,'%% plot forward branch\n'); 
fprintf(fidOut4,'plot(x(end,:),x(1,:),''b-'')\n');
fprintf(fidOut4,'hold on\n\n');
fprintf(fidOut4,'%% plot backward branch\n');
fprintf(fidOut4,'plot(xb(end,:),xb(1,:),''r-'')\n\n');
 
fprintf(fidOut4,'%% plot starting point\n');
fprintf(fidOut4,'plot(xb(end,1),xb(1,1),''k*'')\n');
fprintf(fidOut4,'plot(x(end,1),x(1,1),''r*'')\n');

fclose(fidOut4);


