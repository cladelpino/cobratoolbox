function [model, rxnsInModel, rxnName] = addSinkReactions(model, metabolites, lb, ub)
% Adds a sink reaction for the list of metabolites
%   
% USAGE:
%
%    [model] = addSinkReactions(model, metabolites, lb, ub)
%
% INPUTS:
%    model:          COBRA model structure
%    metabolites:    Cell array of metabolite abreviations as they appear in `model.mets`
%
% OPTIONAL INPUTS:
%    lb:             Lower bounds of reactions
%    ub:             Upper bounds of reactions
%
% OUTPUTS:
%    model:          COBRA model structure containing sink reactions
%    rxnsInModel:    Vector, contains -1 if the reaction did not exist
%                    previously, otherwise it contains the reaction ID of
%                    an identical reaction already present in the model
%
% .. Author: - Ines Thiele 05/06/08
   %         - Claudio Delpino 07/02/19: Modified to use
   %         addMultipleReactions

if ~iscell(metabolites) & ischar(metabolites)
    % assumes it is a single metabolite in a char vector
    metabolites = {metabolites};
end
if any(~ismember(metabolites,model.mets))
    notFound = metabolites(~ismember(metabolites,model.mets)); 
    warning('%s\n','The following metabolites were not found in model and will be added:',notFound{:});
end

rxnsInModel = NaN(size(metabolites));
origMetaboList = metabolites;
% Which metabolites already have a sink reaction ?
% Assumes all single metabolite reactions are sink with a coeff of -1
%rxnIsSink = sum(model.S ~= 0) == 1;
%metHasSink = sum(model.S(:,rxnIsSink),2)<0;
%isLeftOutMet = ismember(metabolites,model.mets) & metHasSink;
%for i = find(isLeftOutMet)'
%   rxnsInModel(strcmp(model.mets{i},metabolites)) = find(model.S(i,:)~=0 & rxnIsSink);
%end

nMets = length(metabolites);
if nargin < 3
    lb = ones(nMets,1)*min(model.lb);
    ub = ones(nMets,1)*max(model.ub);
end

if size(lb,2)==2
    ub = lb(:,2);
    lb = lb(:,1);
end

rxnName = strcat('sink_',metabolites);
k = length(model.rxns) + 1;
model = addMultipleReactions(model,rxnName,metabolites,-eye(nMets),'lb',lb,'ub',ub);
for e = metabolites'
    rxnsInModel(strcmp(e,origMetaboList))=k;
    k = k+1;
end
end
