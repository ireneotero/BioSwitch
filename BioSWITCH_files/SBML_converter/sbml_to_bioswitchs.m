function [] = sbml_to_bioswitchs(sbml_file_name)
% converter from SBML to BioSwitchS format
% requires installation of libSBML-5.15.0-matlab
% example use: sbml_to_bioswitchs('MAPK.xml')
% output: 
%       - MAPK.m -- BioSwitchS input file
%       - MAPK_id_name.txt -- mapping of the species' ids to the
%       compartments they inhabit and their names
% All reactions are assumed of the mass-action type

% Hard-coded name of the zero complex in BioSwitchS. No species should
% have that name. The zero complex is represented as null -> ... or ... -> null in
% the sbml file
complex0NameBioSwitchS = '0';

[pathstr,name,~] = fileparts(sbml_file_name);

if pathstr ~= ''
    pathstr = [pathstr filesep];
end

bioswitchs_file_name = [pathstr name '.m'];
bioswitchs_species_id_name_file_name = [pathstr name '_id_name.txt'];

sbmlobj = TranslateSBML(sbml_file_name); % from libSBML-5.15.0-matlab

% SPECIES -- written in the file
number_of_species = length(sbmlobj.species);
% species = cell(number_of_species,1);
speciesIdName = cell(number_of_species,1);
for i=1:number_of_species
%     species{i} = sbmlobj.species(i).id;
    speciesIdName{i} = [sbmlobj.species(i).id ',' sbmlobj.species(i).compartment ',' sbmlobj.species(i).name];
end

% print a dictionary: species IDs and their corresponding names
fileID = fopen(bioswitchs_species_id_name_file_name,'w');
fprintf(fileID,  '%s',strjoin(speciesIdName,'\n'));
fclose(fileID);


% COMPLEXES
complexes = {};
speciesFromCOM = {};
reactions = cell(length(sbmlobj.reaction),3);
for i=1:length(sbmlobj.reaction)
    reversible = sbmlobj.reaction(i).reversible;
    reactants  = sbmlobj.reaction(i).reactant;
    products   = sbmlobj.reaction(i).product;
    
    % reactants
    currReactants = cell(length(reactants),1);
    for j = 1:length(reactants)
        currReactant = reactants(j).species;
        speciesFromCOM = [speciesFromCOM {currReactant}];
        currStoich   = reactants(j).stoichiometry;
        if currStoich == 1
            currReactants{j} = currReactant;
        else
            currReactants{j} = [int2str(currStoich) '*' currReactant];
        end
    end
    
    currReactants = strjoin(sort(currReactants),'+');
    complexes = [complexes {currReactants}];
    
    % products
    currProducts = cell(length(products),1);
    for j = 1:length(products)
        currProduct = products(j).species;
        speciesFromCOM = [speciesFromCOM {currProduct}];
        currStoich  = products(j).stoichiometry;
        if currStoich == 1
            currProducts{j} = currProduct;
        else
            currProducts{j} = [int2str(currStoich) '*' currProduct];
        end
    end
    
    currProducts = strjoin(sort(currProducts),'+');
    complexes = [complexes {currProducts}];
    
    % save reaction
    reactions{i}{1} = reversible;
    reactions{i}{2} = currReactants;
    reactions{i}{3} = currProducts;
end

species   = unique(speciesFromCOM);
complexes = unique(complexes);
complexes( strcmp(complexes,'') ) = {complex0NameBioSwitchS};


% kinetic constants in BioSwitchS format
kinConsts = {};
zeroComplexId = find(strcmp(complexes,complex0NameBioSwitchS));
for i = 1:length(reactions)
    rev       = reactions{i}{1};
    reactants = reactions{i}{2};
    products  = reactions{i}{3};
    
    if strcmp(reactants,'')
        idReactant = zeroComplexId;
    else
        idReactant = find(strcmp(complexes,reactants));
    end
    
    if strcmp(products,'')
        idProduct = zeroComplexId;
    else
        idProduct = find(strcmp(complexes,products));
    end
    
    currConst = ['k_' int2str(idReactant) '_' int2str(idProduct)];
    kinConsts = [kinConsts {currConst}];
    
    if rev == 1 % reversible reaction -- add another unidirectional raction in the opposite direction
        currConst = ['k_' int2str(idProduct) '_' int2str(idReactant)];
        kinConsts = [kinConsts {currConst}];
    end
end

fileID = fopen(bioswitchs_file_name,'w');
fprintf(fileID, 'function [ESP, COM, KNS] = %s()\n%% Reference: %s\n\n',name,sbml_file_name);
fprintf(fileID, 'ESP={''%s''};\n',strjoin(species,''','''));
fprintf(fileID, 'COM={''%s''};\n',strjoin(complexes,''','''));
fprintf(fileID, 'KNS={''%s''};\n',strjoin(kinConsts,''','''));
fprintf(fileID,'end');
fclose(fileID);
end

