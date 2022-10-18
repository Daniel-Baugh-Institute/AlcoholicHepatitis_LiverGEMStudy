liver_data = readtable('DiseaseData_TPMnorm.txt');
data_struct.genes = liver_data{:,1};
data_struct.tissues = liver_data.Properties.VariableNames(2:end);
data_struct.levels = table2array(liver_data(:,2:end));
data_struct.threshold = 1;

load('Human-GEM.mat')
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
refModel = ihuman;

celltype = [];
hpaData = [];
arrayData = data_struct;
metabolomicsData = [];
removeGenes = true;
taskFile = [];
useScoresForTasks = true;
printReport = true;
taskStructure = essentialTasks;
params = [];
paramsFT = [];

tissuetype=data_struct.tissues;

 for i=1:8
     tissue = char(tissuetype(i)); 
     newGEM = getINITModel2(refModel, tissue, celltype, hpaData,...
         arrayData, metabolomicsData, removeGenes, taskFile, ... 
         useScoresForTasks, printReport, taskStructure, params, paramsFT);
     newGEM.id = char(tissuetype(i)); 
     save(['/GEMGeneration/DiseaseGEM/' char(tissuetype(i)) '.mat'], 'newGEM');
 end

 