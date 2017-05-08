%% Figure out base number of parameters, to then add according to which params will vary across features
% load SAM Structure
% find(cell2mat(arrayfun(@(in1) in1.XSpec.n.n == nPar, SAM.model.variants.tree,'Uni',0)) &&...
%     cell2mat(arrayfun(@(in1) length(in1.XSpec.free.go.freeCat{1}), SAM.model.variants.tree,'Uni',0)))
load('~/perceptualchoice_stop_model/data/2017-05-08/preproc01/subj01/dt5/trialvar/crace/SAM_goTrials')
load('~/perceptualchoice_stop_model/data/2017-03-06/preproc01/subj01/dt5/trialvar/crace/SAM_goTrials')
%%
SAM.model.XCat.name
baseParam = sum(SAM.model.variants.tree(1).XSpec.n.nCat)
SAM.model.variants.tree(1).XSpec.n.nCat
SAM.model.variants.tree(1).XSpec.n.nCatClass

z0Ind = 1;
zCInd = 2;
vInd = 3;
veInd = 4;
t0Ind = 6;
wliwInd = 10;
wffiwInd = 12;

goFeatInd = 1;
stopFeatInd = 2;
stmRow = 1;
rspRow = 2;
cndRow = 3;
%% z0 between responses
nPar = 19;
ind = find(cell2mat(arrayfun(@(in1) in1.XSpec.n.n == nPar, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(z0Ind) == 3, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.features(cndRow,wliwInd,goFeatInd) == 0, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.features(cndRow,wffiwInd,goFeatInd) == 0, SAM.model.variants.tree,'Uni',0)))

%% v and v0 between conditions
nPar = 20;
ind = find(cell2mat(arrayfun(@(in1) in1.XSpec.n.n == nPar, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(vInd) == 3, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(veInd) == 2, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) sum(in1.features(rspRow,:,goFeatInd)) == 0, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.features(cndRow,wliwInd,goFeatInd) == 0, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.features(cndRow,wffiwInd,goFeatInd) == 0, SAM.model.variants.tree,'Uni',0)))

%% z0 between responses, v and v0 between conditions
nPar = 21;
ind = find(cell2mat(arrayfun(@(in1) in1.XSpec.n.n == nPar, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(z0Ind) == 3, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(vInd) == 3, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(veInd) == 2, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) sum(in1.features(rspRow,:,goFeatInd)) == 1, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.features(cndRow,wliwInd,goFeatInd) == 0, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.features(cndRow,wffiwInd,goFeatInd) == 0, SAM.model.variants.tree,'Uni',0)))

%% v and ve between responses, v and v0 between conditions
nPar = 24;
ind = find(cell2mat(arrayfun(@(in1) in1.XSpec.n.n == nPar, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(vInd) == 5, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(veInd) == 4, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.features(cndRow,wliwInd,goFeatInd) == 0, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.features(cndRow,wffiwInd,goFeatInd) == 0, SAM.model.variants.tree,'Uni',0)))

%% z0, v and ve between responses, v and v0 between conditions
nPar = 25;
ind = find(cell2mat(arrayfun(@(in1) in1.XSpec.n.n == nPar, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(z0Ind) == 3, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(vInd) == 5, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.XSpec.n.nCat(veInd) == 4, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.features(cndRow,wliwInd,goFeatInd) == 0, SAM.model.variants.tree,'Uni',0)) & ...
    cell2mat(arrayfun(@(in1) in1.features(cndRow,wffiwInd,goFeatInd) == 0, SAM.model.variants.tree,'Uni',0)))

%%
for i = 1 : length(ind)
% for i = 1 : length(SAM.model.variants.tree)
ind(i)
SAM.model.XCat.name
    SAM.model.variants.tree(ind(i)).XSpec.i.go.iCatClass
    SAM.model.variants.tree(ind(i)).XSpec.i.stop.iCatClass
    SAM.model.variants.tree(ind(i)).XSpec.n.nCat
    SAM.model.variants.tree(ind(i)).XSpec.n.nCatClass
    SAM.model.variants.tree(ind(i)).features

    pause
end
