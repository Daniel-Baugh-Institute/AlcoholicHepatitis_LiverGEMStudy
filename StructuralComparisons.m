% loading in disease-state specific GEM's to perform structural comparison
liver_data = readtable('DiseaseData_TPMnorm.txt');
tissues = liver_data.Properties.VariableNames(2:end);
tissuecat = strcat(tissues,'.mat');

models = {};
for i = 1:8
    load(char(tissuecat(i)));
    models{i} = eval(char(tissues(i)));
end

model_ids = arrayfun(@(i) models{i}.id, (1:numel(models))', 'UniformOutput', false);
res = compareMultipleModels(models);
save('/GEMGeneration/DiseaseGEM/Disease_CompareStruct.mat', 'res');
writecell(res.modelIDs,'/GEMGeneration/DiseaseGEM/DiseaseGEM_modelsIDs.txt')
writecell(res.subsystems.ID,'/GEMGeneration/DiseaseGEM/DiseaseGEM_subsysIDs.txt')
writematrix(res.subsystems.matrix,'/GEMGeneration/DiseaseGEM/DiseaseGEM_subsysmat.txt')
writecell(res.reactions.IDs,'/GEMGeneration/DiseaseGEM/DiseaseGEM_reactionIDs.txt')
writematrix(res.reactions.matrix,'/GEMGeneration/DiseaseGEM/DiseaseGEM_reactionmat.txt')
writematrix(res.structComp,'/GEMGeneration/DiseaseGEM/DiseaseGEM_structComp.txt')