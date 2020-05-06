% The COBRAToolbox: testEfmYieldAnalysis.m
%
% Purpose:
%     - checks yields in a given set of EFMs and (relative) fluxes for desired input and output reactions
%
% Authors:
%     - Created initial test script: Chaitra Sarathy 2 Dec 19
%     - Updates to header docs: Chaitra Sarathy 19 Dec 19

global CBTDIR

% save the current path and initialize the test
currentDir = cd(fileparts(which(mfilename)));

% determine the test path for references
testPath = pwd;

% load the model
model = readCbModel('testModel.mat','modelName','model_irr'); %For all models which are part of this particular test.

% Load reference data
load([testPath filesep 'testData_efmYieldAnalysis.mat']);

% This is a function with three required inputs:
% Case 1: Test whether efmYieldAnalysis gives the expected output when using 3 input and 1 output argument
EFMyield = efmYieldAnalysis(testEFMFluxes, roi_Up, roi_Rel); 
assert(isequal(EFMyield, EFMyield_ref));

% Case 2: Test whether efmYieldAnalysis gives an error when using <3 input and 1 output argument
assert(verifyCobraFunctionError('efmYieldAnalysis', 'inputs', {testEFMFluxes}, 'outputArgCount', 1));
assert(verifyCobraFunctionError('efmYieldAnalysis', 'inputs', {testEFMFluxes, roi_Up}, 'outputArgCount', 1));

