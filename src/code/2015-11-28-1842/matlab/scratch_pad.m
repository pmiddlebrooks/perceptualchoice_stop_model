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


%%
% load SAM Structure
% find(cell2mat(arrayfun(@(in1) in1.XSpec.n.n == nPar, SAM.model.variants.tree,'Uni',0)) &&...
%     cell2mat(arrayfun(@(in1) length(in1.XSpec.free.go.freeCat{1}), SAM.model.variants.tree,'Uni',0)))
load('SAM_goTrials')
baseParam = sum(SAM.model.variants.tree(1).XSpec.n.nCat)
SAM.model.variants.tree(1).XSpec.n.nCat
SAM.model.variants.tree(1).XSpec.n.nCatClass

%%

nPar = 28;

ind = find(cell2mat(arrayfun(@(in1) in1.XSpec.n.n == nPar, SAM.model.variants.tree,'Uni',0)))
% ind = find(cell2mat(arrayfun(@(in1) sum(in1.XSpec.n.nCat == nPar), SAM.model.variants.tree,'Uni',0)))
% SAM.model.variants.tree(9).XSpec.n.nCat

%%
for i = 1 : length(ind)
% for i = 1 : length(SAM.model.variants.tree)
ind(i)
    SAM.model.variants.tree(ind(i)).XSpec.i.go.iCatClass
    SAM.model.variants.tree(ind(i)).XSpec.n.nCat
    SAM.model.variants.tree(ind(i)).XSpec.n.nCatClass

    pause
end
    