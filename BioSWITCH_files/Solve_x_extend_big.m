function [k_extend, x_extend, set_idx, Str_sol] = Solve_x_extend_big(ESP,COM,KNS,network_file)

 %warning('off','symbolic:sym:isAlways:TruthUnknown')      
         
out=Deficiency_computations(ESP,COM,KNS);


lhs=out.A_kPsi;
rhs=out.rhs;


Psi=out.Psi; 
Psi_c=out.Psi_c;

alpha_vec=out.alpha_vec;

ESP=out.esp;
COM=out.com;
KNS=out.kns;

n_com=out.n_com;
n_link=out.n_link;
n_esp=out.n_esp;
n_rec=out.n_rec;

def=out.def;

matrix_N=matrixN(ESP,COM,KNS);

s=rank(matrix_N);

lambda=n_esp-s;

% Compute H and Hs
[reduced_lhs, reduced_rhs]=Extract_Lin_Indep_Eqs(lhs,rhs,Psi);

for ii=1:1:n_com-n_link
   H(ii)=reduced_lhs(ii)-reduced_rhs(ii);
   H_s(ii)=subs(H(ii),Psi,Psi_c);
end

if size(H,2)>size(H,1)
    H=H.';
end

if size(H_s,2)>size(H_s,1)
    H_s=H_s.';
end


for ii=1:1:n_esp
    eval(sprintf('syms x_%d',ii));
    eval(sprintf('c(%d)=x_%d;',ii,ii));
end

if size(c,2)>size(c,1)
    c=c.';
end


% Decision variables
x_extend(1:n_esp)=c;
x_extend(n_esp+1:n_esp+def)=alpha_vec;


nk=size(H_s,1);

name_file = sprintf('%s_manifold.txt',network_file);
    fidOut=fopen(name_file,'w');
    
    fprintf(fidOut,'Manifold equations\n\n');
    
    for ii=1:1:n_com-n_link
        fprintf(fidOut,'%s=0\n',char(H_s(ii)));
    end
    
    fprintf(fidOut,'\n');
    fclose(fidOut)
    
% if n_esp >= nk  % Proper and under-dimensioned    
%     ncomb = factorial(n_esp)/factorial(n_esp-nk)/factorial(nk); %IOM check
%     indices = 1:1:n_esp;       
% else % over-dimensioned   
%     ncomb = factorial(n_esp+def)/factorial(n_esp+def-nk)/factorial(nk); %IOM check
%     indices = 1:1:n_esp+def;     
% end
    
    indices2 = 1:1:def;
    nk2 = def + n_esp - nk;
    
    
    set2 = nchoosek(indices2, nk2);
    ncomb2 = factorial(def)/factorial(def-nk2)/factorial(nk2);
    
   % set = nchoosek(indices,nk);
   
   
    ncomb = round(ncomb2);
    Str_sol={};
    set_idx=[];
    
    ii=1;
    
%     name_file = sprintf('%s_manifold.txt',network_file);
%     fidOut=fopen(name_file,'w');
%     
%     fprintf(fidOut,'Manifold equations\n\n');
%     
%     for ii=1:1:n_com-n_link
%         fprintf(fidOut,'%s=0\n',char(H_s(ii)));
%     end
%     
%     fprintf(fidOut,'\n');
%     
%       fidOut=fopen(name_file,'w');
%     fprintf(fidOut,'Var Indep\n\n');
%     
%       
%     for ii=1:1:ncomb
%         fprintf(fidOut,'var_indep_%d=%s\n',ii,char(x_extend(set(ii,:))));
%     end
%     
%     fclose(fidOut);
    
    %
    
    
    sprintf('Selection of independent variables: Choose 0 for automatic selection, Choose 1 for manual selection')
    selec = input('Type 0 (automatic), 1 (manual):');
    if selec == 0
        ii=299;
        while ii<=ncomb2 && isempty(Str_sol)
            ii
            vars=setxor(x_extend,alpha_vec(set2(ii,:))).';
            vars
          
            S = solve(H_s==0, vars);
            S
            
            EmptyIndex = structfun(@isempty,S);
            idi=find(EmptyIndex>0);
            if isempty(idi)
                Str_sol=S;
                set_idx=set2(ii,:); %OJO
            end
            ii=ii+1;
            if ii>=ncomb2 && isempty(Str_sol)
                error('Matlab solve was not able to find the expression of the manifold equations (stored in txt file), you still can try by hand')
                break
            end
        end
        
    elseif selec == 1
%         sprintf('There are %d possible vectors of independent variables (check %s_manifold.txt), choose a vector from 1 to %d', ncomb, network_file, ncomb )
%         selectopt = input('chosen vector:');
%         %
%         x_extend(set(selectopt,:));
%         S{selectopt} = solve(H_s==0,x_extend(set(selectopt,:)));
%         S{selectopt}
%         %
%         EmptyIndex = structfun(@isempty,S{selectopt});
%         idi=find(EmptyIndex>0);
%         if isempty(idi)
%             Str_sol=S{selectopt};
%             set_idx=set(selectopt,:);
%         else
%             error('Matlab solve was not able to find the expression of the manifold equations (stored in txt file) for the vector chosen, try another or use the automatic selection')
%             return
%         end
%         %     end
        %
       
    end
    
    
    %  When ncomb  = 1, solve(H_s==0,zeros(nk,1),x_extend(set(ii,:))); gives
    % the error "The number of rows of the string matrix must match the
    % number of elements in the cell." and does not give a symbolic
    % solution.
    
    
    
    for ii=1:1:n_rec
        eval(sprintf('syms %s',KNS{ii}));
        eval(sprintf('k(%d)=%s;',ii,KNS{ii}));
    end
    
    if size(k,2)>size(k,1)
        k=k.';
    end
    
    k_extend(1:n_rec)=k;
    k_extend(n_rec+1:n_rec+lambda)=setxor(x_extend,x_extend(set_idx));
    
    if size(k_extend,2)>size(k_extend,1)
        k_extend=k_extend.';
    end

    
