%% SPOT and Eflux2  %%

addpath(genpath('GEMGeneration'))
addpath(genpath('FluxBalanceAnalysis'))

%% For Disease-state specific GEMS %%

% reading in the data
liver_data = readtable('DiseaseData_TPMnorm.txt');
tissues = liver_data.Properties.VariableNames(2:end);
tissuecat = strcat(tissues,'.mat');

models = {};
for i = 1:8
    load(char(tissuecat(i)));
    models{i} = eval(char(tissues(i)));
end


obj = {};
flux = {};
model = {};
for i=1:8
    model = char(tissuecat(i));
    [obj_func, eflux2_flux, fbamodel] = SpotEflux2(i,model,'DiseaseData_TPMnorm.txt');
    obj{i} = obj_func;
    flux{i} = eflux2_flux;
    model{i} = fbamodel;
end

save('FluxBalanceAnalysis/DiseaseState_objfunc.mat','obj');
save('FluxBalanceAnalysis/DiseaseState_flux.mat','flux');
save('FluxBalanceAnalysis/DiseaseState_models.mat','model');

ihuman_rxns = readtable('Human-GEM-reactions.txt','Delimiter','\t');
ihuman_rxns = table2array(ihuman_rxns);
load('Human-GEM.mat') %saves as ihuman
ihuman_subsys = ihuman.subSystems;

allfluxes = zeros(size(ihuman_rxns,1),8);
% putting all fluxes in a matrix and saving to text file
for i =1:8
    for j = 1:size(ihuman_rxns,1)
        findrxn = find(strcmp(char(ihuman_rxns(j)),model{i}.rxns));
        if isempty(findrxn)
            findflux = 0;
        else
            findflux = flux{i}(findrxn);
        end
        allfluxes(j,i) = findflux;
    end  
end

allfluxes = array2table(allfluxes);
allfluxes = [ihuman_rxns,ihuman_subsys, allfluxes];
allfluxes.Properties.VariableNames = {'Reaction','Subsystem', tissues}; 

allfluxes_gt01 = [];
for i = 1:size(ihuman_rxns,1)
    if any(abs(table2array(allfluxes(i,3:10)))>0.1)
        allfluxes_gt01 = [allfluxes(i,:); allfluxes_gt01];
    end
end

allfluxes_gt01.Properties.VariableNames = {'Reaction','Subsystem', tissues}; 

writetable(allfluxes,'FluxBalanceAnalysis/DiseaseStateGEM_allfluxes.txt','Delimiter','\t')
writetable(allfluxes_gt01,'FluxBalanceAnalysis/DiseaseStateGEM_fluxes_gt0.1.txt','Delimiter','\t')

