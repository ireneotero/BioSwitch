function [Yr, Nr] = matrixYrNr(ESP,COM,KNS,input_species_vector)
%[Yr, Nr] = matrixYrNr(ESP,COM,KNS,input_species)

n_esp = size(ESP,2);

n_rec = size(KNS,2);

n_rec_inner = n_rec;

input_species = input_species_vector;

for ii = 1:1:n_esp;
    species.name{ii}=ESP{ii};
    species.index(ii)=ii;
end


for ii = 1:1:n_rec;
    reaction.kname{ii}=KNS{ii};
    reaction.index(ii)=ii;
end


Y = matrixY(ESP,COM);

% check species-alone complexes in the inner network

species_alone_counter = 0;

for ii=1:1:size(Y,2)
   sum_col(ii) = sum(Y(:,ii));
   if sum_col(ii)==1;
       species_alone_counter = species_alone_counter+1;
       species_alone_complexes.name{species_alone_counter}=COM{ii};
       species_alone_complexes.complex_index(species_alone_counter)=ii;      
       species_alone_complexes.species_index(species_alone_counter)=find(strcmp(species_alone_complexes.name{species_alone_counter},species.name));
   end       
end

% add complex 0

index_complex_0 = size(COM,2)+1;

COM{index_complex_0} = '0';


% add degradation reactions for not-alone species
species_not_alone.species_index = setxor(1:1:n_esp,species_alone_complexes.species_index);


for ii=1:1:max(size(species_not_alone.species_index))
   eval(sprintf('COM{index_complex_0+ii}=''%s'';',ESP{species_not_alone.species_index(ii)}));
   species_not_alone.complex_index(ii)=index_complex_0+ii;
end

counter1=0;
for ii=1:1:max(size(species_not_alone.species_index))
    counter1=counter1+1;
   reaction.index(n_rec+counter1)=n_rec+counter1;
   reaction_kname=sprintf('k_%d_%d',species_not_alone.complex_index(ii),index_complex_0);
   reaction.kname{n_rec+counter1}=reaction_kname;
end
n_rec=n_rec+counter1;

% add degradation reactions for species-alone complexes in the inner network
counter2=0;
for ii=1:1:size(species_alone_complexes.species_index,2)
    counter2=counter2+1;
    reaction.index(n_rec+counter2)=n_rec+counter2;
    reaction_kname = sprintf('k_%d_%d',species_alone_complexes.complex_index(ii),index_complex_0);
    reaction.kname{n_rec+counter2}=reaction_kname;
end
n_rec=n_rec+counter2;
n_rec_output=counter1+counter2;


% add inflow reactions for input species-alone
counter3=0;
for ii=1:1:max(size(input_species))
    for jj=1:1:size(species_alone_complexes.species_index,2)
      iden(ii,jj) = isequal(input_species(ii), species_alone_complexes.species_index(jj));
      if iden(ii,jj)>0
          counter3=counter3+1;
          reaction.index(n_rec+counter3)=n_rec+counter3;
          reaction_kname = sprintf('k_%d_%d',index_complex_0,species_alone_complexes.complex_index(jj));
          reaction.kname{n_rec+counter3}=reaction_kname;
      end
    end
end
n_rec=n_rec+counter3;

% add inflow reactions for input species-not-alone
counter4=0;
for ii=1:1:max(size(species_not_alone.species_index))
 counter4=counter4+1;
 reaction.index(n_rec+counter4)=n_rec+counter4;
 reaction_kname=sprintf('k_%d_%d',index_complex_0,species_not_alone.complex_index(ii));
 reaction.kname{n_rec+counter4}=reaction_kname;
end
n_rec=n_rec+counter4;
n_rec_input=counter3+counter4;


for ii=size(KNS,2)+1:n_rec
     eval(sprintf('KNS{ii}=''%s'';',reaction.kname{ii}));
end


% Find source complex for reactions
for ii=1:1:n_rec_inner   
    cont=3;       
        while reaction.kname{ii}(cont)~='_'         
           cte_so{ii}(cont-2)=reaction.kname{ii}(cont);
           cont=cont+1;
        end
       
        eval(sprintf('reaction.source_index(%d)=%d;',ii,str2num(cte_so{ii})))
       
end


% Find sink complex for reactions
   
for ii=1:1:n_rec_inner    
    cont=max(strfind(reaction.kname{ii},'_'))+1; 
            while cont<=size(reaction.kname{ii},2)
           cte_si{ii}(cont-max(strfind(reaction.kname{ii},'_')))=reaction.kname{ii}(cont);
           cont=cont+1;       
          end        
        eval(sprintf('reaction.sink_index(%d)=%d;',ii,str2num(cte_si{ii})));
end




% Compute matrix Yr and Nr
for ii=1:1:n_rec_inner  
 Yr(:,ii)=Y(:,reaction.source_index(ii));
end

for ii=1:1:n_rec_inner
 Nr(:,ii)=Y(:,reaction.sink_index(ii))-Y(:,reaction.source_index(ii));
end
 
Yr=[Yr eye(n_esp)];
Nr=[Nr -eye(n_esp)];






