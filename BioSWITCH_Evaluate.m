function BioSWITCH_Evaluate(network_file)
% BioSWITCH_Evaluate(network_file)
% evaluates the network (rank of stoichometric matrix, conservation laws,
% deficiency)

if nargin < 1
    error('Name of network file (char) needs to be provided as an argument')
end

%eval(sprintf('[ESP, COM, KNS] = %s;',network_file));
network_file_function = str2func(network_file);
[ESP, COM, KNS] = network_file_function();

n_esp=size(ESP,2);

n_reac=size(KNS,2);

n_com=size(COM,2);


matrix_N = matrixN(ESP,COM,KNS);
n_link = nlink(ESP,COM,KNS);
B = Compute_B(ESP,COM,KNS);


rankN = rank(matrix_N);


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

if size(b,2)>size(b,1)
    b=b.';end
end

n_cons_laws = n_esp - double(rankN);

nfile = sprintf('Evaluation_%s.txt',network_file);
FID=fopen(nfile,'w+');


fprintf('\n')
fprintf(FID,'The number of species (n) is: %d\n\n', n_esp);
fprintf('The number of species (n) is: %d\n\n', n_esp);

for ii=1:1:n_esp
    fprintf(FID,'species %d: %s \n', ii, ESP{ii});
    fprintf('species %d: %s \n', ii, ESP{ii});
end

fprintf(FID,'\n');
fprintf('\n');
fprintf(FID,'The number of complexes (m) is: %d\n\n', n_com);
fprintf('The number of complexes (m) is: %d\n\n', n_com);

for ii=1:1:n_com
    fprintf(FID,'complex %d: %s \n', ii, COM{ii});
    fprintf('complex %d: %s \n', ii, COM{ii});
end

fprintf(FID,'\nThe number of linkage (l) classes is: %d\n',double(n_link));
fprintf('\nThe number of linkage classes (l) is: %d\n',double(n_link));

if n_esp == n_com-n_link
    fprintf(FID,'\nn = m-l (proper network)\n');
    fprintf('\nn = m-l (proper network)\n');
elseif n_esp > n_com-n_link
    fprintf(FID,'\nn > m-l (underdimensioned network)\n');
    fprintf('\nn > m-l (underdimensioned network)\n');
elseif n_esp < n_com-n_link
    fprintf(FID,'\nn < m-l (overdimensioned network)\n');
    fprintf('\nn < m-l (overdimensioned network)\n');
end

fprintf(FID,'\nThe rank of the stoichiometrix matrix (rho) is: %d\n',double(rankN));
fprintf('\nThe rank of the stoichiometrix matrix (rho) is: %d\n',double(rankN));
fprintf(FID,'\nThe number of mass conservation laws (n-rho) is: %d\n', n_cons_laws);
fprintf('\nThe number of mass conservation laws (n-rho) is: %d\n', n_cons_laws);

if n_cons_laws ~=  0
    fprintf(FID,'\n');
    fprintf('\n');
    for ii=1:1:size(b,1)
    fprintf(FID,'Conservation law %d: %s\n',ii,char(b(ii)));
    fprintf('Conservation law %d: %s\n',ii,char(b(ii)));
    end
else
    fprintf(FID,'\n');
    fprintf('\n');
    fprintf(FID,'The network has no mass conservation laws');
    fprintf('The network has no mass conservation laws');
end


id_complex_0 = find(strcmp(COM,'0'));


for ii=1:1:n_reac  
    k.name{ii} = KNS{ii};
end
 
for ii=1:1:n_reac
    str1 = strcat('_',num2str(id_complex_0),'_');
    str2 = strcat('_',num2str(id_complex_0));
    pattend = (length(KNS{ii}) >= length(str2) && strcmp(KNS{ii}(length(KNS{ii})-length(str2)+1:length(KNS{ii})), str2));
    
    if  ~isempty(strfind(KNS{ii}, str1))
        k.type{ii} = 'inflow reaction';
    elseif  isempty(strfind(KNS{ii}, str1)) && pattend
        k.type{ii} ='outflow reaction';
    else
        k.type{ii} = 'true reaction';
    end
end

fprintf(FID,'\n');
fprintf('\n');
fprintf(FID,'Reactions are classified in:\n');
fprintf('Reactions are classified in:\n');

for ii=1:1:n_reac
    fprintf(FID,'%s : %s\n',k.name{ii},k.type{ii});
    fprintf('%s : %s\n',k.name{ii},k.type{ii});
end

outd=Deficiency_computations(ESP,COM,KNS)
deficiency = outd.def;

fprintf(FID,'\n');
fprintf('\n');
fprintf(FID,'The deficiency of the network is: %d\n',deficiency);
fprintf('The deficiency of the network is: %d\n',deficiency);

if deficiency == 0
    fprintf(FID,'If the network is weakly reversible: there is a unique steady state per stoichiometric compatibility class (for all parameter sets)\n');
    fprintf(FID,'If the network is not weakly reversible, no positive equilibrium is admitted \n');
    fprintf('If the network is weakly reversible: there is a unique steady state per stoichiometric compatibility class (for all parameter sets)\n');
    fprintf('If the network is not weakly reversible, no positive equilibrium is admitted \n');
end

 fprintf('\n');
 
fclose(FID);  
