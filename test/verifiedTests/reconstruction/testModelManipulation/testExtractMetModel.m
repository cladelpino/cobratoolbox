% The COBRAToolbox: testUpdateGenes.m
%
% Purpose:
%     - testExtractMetModel tests extractMetModel
%
% Authors:
%     - Uri David Akavia August 2017

global CBTDIR

% save the current path
currentDir = pwd;

fileDir = fileparts(which('testUpdateGenes'));
cd(fileDir);

if 1
    model = getDistributedModel('Recon2.v05.mat');
else
    if 0
        model = getDistributedModel('Recon2.v04.mat');
    else
        load('Recon2.v04.mat');
        model = modelR204;
        %delete gpr for ATPS4m
        model.grRules(strcmp(model.rxns,'ATPS4m')) = {''};
        model.rules(strcmp(model.rxns,'ATPS4m')) = {''};
        model = convertOldStyleModel(model, 0);
        res = verifyModel(model, 'silentCheck', true);
    end
end
                
%compare against explicitly loaded models to conserve the ids.
load('testExtractMetModel.mat', 'emptyModel', 'atpModel', 'pppg9Level0', 'pppg9Level1');

% Test getting level 0 (just reactions that involve a metabolite)
model2 = extractMetModel(model, 'pppg9', 0, 1);

printLevel=1;
pppg9Level0.modelID=model.modelID;
[isSame, nDiff, commonFields] = isSameCobraModel(model2, pppg9Level0, printLevel);
assert(isSame);

% Test getting level 1 (include one reaction away from reactions that involve a metaoblite)
model2 = extractMetModel(model, 'pppg9', 1, 1);
pppg9Level1.modelID=model.modelID;
assert(isSameCobraModel(model2, pppg9Level1));

% Test asking for a very common metabolite, empty model should be returned
model2 = extractMetModel(model, 'atp', 0, 1);
emptyModel.modelID=model.modelID;
[isSame, nDiff, commonFields] = isSameCobraModel(model2, emptyModel, printLevel);
assert(isSame);

% Test asking for a very common metabolite, with high limit on connectivity
model2 = extractMetModel(model, 'atp', 0, 1, 99999);
atpModel.modelID=model.modelID;
atpModel.grRules=model2.grRules;
atpModel.rules=model2.rules;
[isSame, nDiff, commonFields] = isSameCobraModel(model2, atpModel, printLevel);
assert(isSame);

%return to original directory
cd(currentDir)