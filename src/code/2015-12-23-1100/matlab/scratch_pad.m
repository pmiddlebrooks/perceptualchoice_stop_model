%%
load('data/raw/broca_behavior2.mat')

%%
ssdList         = unique(trialData.ssd(~isnan(trialData.ssd), :));
nSSDMin         = 100; % minimum # of stop trials for a SSD to be included

nSSD            = nan(length(ssdList), 1);

for i = 1 : length(ssdList)
    nSSD(i) = sum(trialData.ssd == ssdList(i));
end

removeSSD = ssdList(nSSD < nSSDMin);
Lia = ismember(trialData.ssd, removeSSD);
trialData = trialData(~Lia, :);    
%%
io.rootDir                    = fullfile(fileparts(fileStr.preprocDataDir),'subj01');
io.behavFile                  = fullfile(io.rootDir,'data_subj01.mat');
load(io.behavFile)


%% Figure out base number of parameters, to then add according to which params will vary across features
% load SAM Structure
% find(cell2mat(arrayfun(@(in1) in1.XSpec.n.n == nPar, SAM.model.variants.tree,'Uni',0)) &&...
%     cell2mat(arrayfun(@(in1) length(in1.XSpec.free.go.freeCat{1}), SAM.model.variants.tree,'Uni',0)))
load('SAM_goTrials')
baseParam = sum(SAM.model.variants.tree(1).XSpec.n.nCat)
SAM.model.variants.tree(1).XSpec.n.nCat
SAM.model.variants.tree(1).XSpec.n.nCatClass

%%
z0Ind = 1;
zCInd = 2;
vInd = 3;
veInd = 4;
t0Ind = 6;

nPar = 18 % baseline numbe of parameters
ind = find(cell2mat(arrayfun(@(in1) in1.XSpec.n.n == nPar, SAM.model.variants.tree,'Uni',0)))
% ind = find(cell2mat(arrayfun(@(in1) sum(in1.XSpec.n.nCat == nPar), SAM.model.variants.tree,'Uni',0)))
% SAM.model.variants.tree(9).XSpec.n.nCat


%% z0 between responses, v and v0 between conditions
nPar = 23;
ind = find(cell2mat(arrayfun(@(in1) in1.XSpec.n.n == nPar, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(z0Ind) == 3, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(vInd) == 4, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(veInd) == 3, SAM.model.variants.tree,'Uni',0)))

%% z0 and t0 between responses, v and v0 between conditions
nPar = 24;
ind = find(cell2mat(arrayfun(@(in1) in1.XSpec.n.n == nPar, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(z0Ind) == 3, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(t0Ind) == 3, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(vInd) == 4, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(veInd) == 3, SAM.model.variants.tree,'Uni',0)))


%% z0 and zC between responses, v and v0 between conditions
nPar = 24;
ind = find(cell2mat(arrayfun(@(in1) in1.XSpec.n.n == nPar, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(z0Ind) == 3, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(zCInd) == 3, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(vInd) == 4, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(veInd) == 3, SAM.model.variants.tree,'Uni',0)))


%%
for i = 1 : length(ind)
% for i = 1 : length(SAM.model.variants.tree)
ind(i)
    SAM.model.variants.tree(ind(i)).XSpec.i.go.iCatClass
    SAM.model.variants.tree(ind(i)).XSpec.n.nCat
    SAM.model.variants.tree(ind(i)).XSpec.n.nCatClass
    SAM.model.variants.tree(ind(i)).features

    pause
end
    