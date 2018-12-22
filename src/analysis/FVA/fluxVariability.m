function [minFlux, maxFlux, Vmin, Vmax] = fluxVariability(model, optPercentage, osenseStr, rxnNameList, printLevel, allowLoops, method, cpxControl, advind)
% Performs flux variablity analysis
%
% USAGE:
%
%    [minFlux, maxFlux] = fluxVariability(model, optPercentage, osenseStr, rxnNameList, printLevel, allowLoops, method)
%
% INPUT:
%    model:            COBRA model structure
%
% OPTIONAL INPUTS:
%    optPercentage:    Only consider solutions that give you at least a certain
%                      percentage of the optimal solution (Default = 100
%                      or optimal solutions only)
%    osenseStr:        Objective sense ('min' or 'max') (Default = 'max')
%    rxnNameList:      List of reactions for which FVA is performed
%                      (Default = all reactions in the model)
%    printLevel:       Verbose level (default: 0)
%    allowLoops:       Whether loops are allowed in solution. (Default = true)
%                      See `optimizeCbModel` for description
%                      Choose different values for allowLoops to use different methods for loopless FVA (implemented in Nov 2017):
%                           *  0 : the original loopless FVA
%                           * -1 : loopless FVA with Fast-SNP preprocessing of nullspace
%                           * -2 : localized loopless FVA using information from nullsapce
%                           * -3 : localized loopless FVA using information from EFMs. 
%                                  Require CalculateFluxModes.m from EFMtool to calculate EFMs.
%    method:           when Vmin and Vmax are in the output, the flux vector can be (Default = 2-norm):
%
%                        * 'FBA'    : standards FBA solution
%                        * '0-norm' : minimzes the vector  0-norm
%                        * '1-norm' : minimizes the vector 1-norm
%                        * '2-norm' : minimizes the vector 2-norm
%                        * 'minOrigSol' : minimizes the euclidean distance of each vector to the original solution vector
%
%   cpxControl:        solver-specific parameter structure
%
%   advind:            switch to use the solution basis
%
%                           - 0 : default
%                           - 1 : uses the original problem solution basis as advanced basis
%
% OUTPUTS:
%    minFlux:          Minimum flux for each reaction
%    maxFlux:          Maximum flux for each reaction
%
% OPTIONAL OUTPUT:
%    Vmin:             Matrix of column flux vectors, where each column is a
%                      separate minimization.
%    Vmax:             Matrix of column flux vectors, where each column is a
%                      separate maximization.
%
% .. Authors:
%       - Markus Herrgard  8/21/06 Original code.
%       - Ronan Fleming   01/20/10 Take the extremal flux from the flux vector,
%                         not from the objective since this is invariant
%                         to the value and sign of the coefficient
%       - Ronan Fleming   27/09/10 Vmin, Vmax
%       - Marouen Ben Guebila 22/02/2017 Vmin,Vmax method

global CBT_LP_PARAMS

if nargin < 2 || isempty(optPercentage)
    optPercentage = 100;
end
if nargin < 3 || isempty(osenseStr)
    if isfield(model, 'osenseStr')
        osenseStr = model.osenseStr;
    elseif isfield(model, 'osense')
        if model.osense == 1
            osenseStr = 'min';
        elseif model.osense == -1
            osenseStr = 'max';
        else
            error('model.osense can only be either 1 (for minimization) or -1 (for maximization)')
        end
    else
        % default max
        osenseStr = 'max';
    end
end
if nargin < 4 || isempty(rxnNameList)
    rxnNameList = model.rxns;
end
if nargin < 5
    printLevel = 0;
end
if nargin < 6
    allowLoops = true;
end
if nargin < 7
    method = '2-norm';
end
if nargin < 8
    cpxControl = struct();
end
if nargin < 9
   advind = 0;
end

%Stop if there are reactions, which are not part of the model
if any(~ismember(rxnNameList,model.rxns))
    presence = ismember(rxnNameList,model.rxns);
    error('There were reactions in the rxnList which are not part of the model:\n%s\n',strjoin(rxnNameList(~presence),'\n'));
end

% Set up the problem size
[nMets,nRxns] = size(model.S);
Vmin=[];
Vmax=[];
if nargout > 2
    OutputMatrix = 1;
else
    OutputMatrix = 0;
end

% LP solution tolerance
if exist('CBT_LP_PARAMS', 'var')
    if isfield(CBT_LP_PARAMS, 'objTol')
        tol = CBT_LP_PARAMS.objTol;
    else
        tol = 1e-6;
    end
    if nargout < 3
        minNorm = 0;
    else
        minNorm = 1;
    end
end

%Return if minNorm is not FBA but allowloops is set to false
%This is currently not supported as it requires mechanisms that are likely
%incompatible.
if allowLoops <= 0 && minNorm && ~strcmp(method,'FBA')
    error('Cannot return solutions with special properties if allowLoops is set to false.\nIf you want solutions without loops please set method to ''FBA''.');
end
% Determine constraints for the correct space (0-100% of the full space)
if sum(model.c ~= 0) > 0
    hasObjective = true;
else
    hasObjective = false;
end

if printLevel == 1
    showprogress(0,'Flux variability analysis in progress ...');
end
if printLevel > 1
    fprintf('%4s\t%4s\t%10s\t%9s\t%9s\n','No','Perc','Name','Min','Max');
end

if ~isfield(model,'b')
    model.b = zeros(size(model.S,1),1);
end
% Set up the general problem
rxnListFull = model.rxns;
LPproblem.c = model.c;
LPproblem.lb = model.lb;
LPproblem.ub = model.ub;
if ~isfield(model,'csense')
    LPproblem.csense(1:nMets,1) = 'E';
else
    LPproblem.csense = model.csense(1:nMets);

    % print a warning message if the csense vector does not have the same length as the mets vector
    if length(model.mets) ~= length(model.csense)
        warning(' > model.csense does not have the same length as model.mets. Consider checking the model using >> verifyModel.');
    end
end
LPproblem.csense = columnVector(LPproblem.csense);
LPproblem.A = model.S;
LPproblem.b = model.b;

%solve to get the original model optimal objective
if strcmp(osenseStr, 'max')
    LPproblem.osense = -1;
else
    LPproblem.osense = 1;
end

% Solve initial (normal) LP
loopInfo = struct();
if allowLoops == 1
    tempSolution = solveCobraLP(LPproblem, cpxControl);
else
    % Both Fast-SNP (-1) and solving an MILP (-2) return a minimal feasible nullspace
    if allowLoops == 0
        % find the usual internal nullspace (Schellenberger et al., 2009)
        [MILPproblem, loopInfo] = addLoopLawConstraints(LPproblem, model, 1:nRxns, 1);
    elseif allowLoops == -1
        % find a minimal feasible nullspace Fast-SNP (Saa and Nielson, 2016)
        if printLevel
            fprintf('Reduce complexity by nullspace preprocessing (Fast-SNP)\n')
        end
        [MILPproblem, loopInfo] = addLoopLawConstraints(LPproblem, model, 1:nRxns, 2);
    elseif allowLoops <= -2
        % find a minimal feasible nullspace by one single MILP (Chan et al., 2017)
        % AND implement localized loopless constraints
        if printLevel
            fprintf('Reduce complexity by nullspace preprocessing and implementing localized loopless constraints (LLCs)\n')
        end
        useRxnLink = false;
        if allowLoops == -2
            [MILPproblem, loopInfo] = addLoopLawConstraints(LPproblem, model, 1:nRxns, 3);
        elseif allowLoops <= -3
            [MILPproblem, loopInfo] = addLoopLawConstraints(LPproblem, model, 1:nRxns, 4);
            if isempty(loopInfo.rxnLink)
                if printLevel
                    fprintf('Unable to find EFMs. Use connections from nullspace to implement LLCs\n')
                end
            else
                useRxnLink = true;
                if printLevel
                    fprintf('Use connections from EFMs to implement LLCs\n')
                end
            end
        end
        [alwaysLLC, rxnInLoopsAlwaysOn, conCompAlwaysOn, x0] = preprocessLLC(LPproblem, ...
            model, nRxns, loopInfo.rxnInLoops, osenseStr, loopInfo.conComp, printLevel);        
        if alwaysLLC
            % apply LLCs for loopless FBA if the objective function or
            % additional constraints contain reaction fluxes that require
            % loop law for optimality under the loopless condition
            MILPproblem = updateLLCs(MILPproblem, conCompAlwaysOn, rxnInLoopsAlwaysOn, loopInfo, [], useRxnLink);
        end
    end
    tempSolution = solveCobraMILP(MILPproblem);
end

if tempSolution.stat == 1
    if strcmp(osenseStr,'max')
        objValue = floor(tempSolution.obj/tol)*tol*optPercentage/100;
    else
        objValue = ceil(tempSolution.obj/tol)*tol*optPercentage/100;
    end
else
    error('The FVA could not be run because the model is infeasible or unbounded')
end

%set the objective
if hasObjective
    LPproblem.A = [model.S;columnVector(model.c)'];
    LPproblem.b = [model.b;objValue];
    if strcmp(osenseStr, 'max')
        LPproblem.csense(end+1) = 'G';
    else
        LPproblem.csense(end+1) = 'L';
    end
end

%get the initial basis
if advind == 1
    LPproblem.basis = tempSolution.basis;
end
LPproblem.S = LPproblem.A;%needed for sparse optimisation

% Loop through reactions
maxFlux = zeros(length(rxnNameList), 1);
minFlux = zeros(length(rxnNameList), 1);

%Thats not true. The Euclidean norm does not get rid of loops if the
%objective reaction is part of the loop.
% if length(minNorm)> 1 || minNorm > 0
%     %minimizing the Euclidean norm gets rid of the loops, so there
%     %is no need for a second slower MILP approach
%     allowLoops=1;
% end

solutionPool = zeros(length(model.lb), 0);

v=ver;
PCT = 'Parallel Computing Toolbox';
if  any(strcmp(PCT,{v.Name})) && license('test','Distrib_Computing_Toolbox')
    p = gcp('nocreate');
    if isempty(p)
        poolsize = 0;
    else
        poolsize = p.NumWorkers;
    end
    PCT_status=1;
else
    PCT_status=0;  % Parallel Computing Toolbox not found.
end

minFlux = model.lb(ismember(model.rxns,rxnNameList));
maxFlux = model.ub(ismember(model.rxns,rxnNameList));
preCompMaxSols = cell(nRxns,1);
preCompMinSols = cell(nRxns,1);

%We will calculate a min and max sum flux solution.
%This solution will (hopefully) provide multiple solutions for individual
%reactions.
QuickProblem = LPproblem;
[Presence,Order] = ismember(rxnNameList,model.rxns);
QuickProblem.c(:) = 0;
QuickProblem.c(Order(Presence)) = 1;
if allowLoops <= 0 && allowLoops >= -1
    % map allowLoops = 0, -1 to preprocessMethod = 1 (original), 2 (Fast-SNP)
    % Skip this when using localized loopless constraints (LLCs) for two reasons:
    % 1. With loopless constraints, one does not expect to have many reactions hitting the bounds
    % 2. LLCs invoke a subset of all loopless constraints and the associated binary
    %    variables. Maximize everything at once will require invoking almost all binary variables
    QuickProblem = addLoopLawConstraints(QuickProblem, model, 1:nRxns, -allowLoops + 1, loopInfo);
end
%Maximise all reactions
QuickProblem.osense = -1;
if allowLoops <= 0
    if allowLoops >= -1
        sol = solveCobraMILP(QuickProblem);
    else
        sol.full = NaN(nRxns, 1);
    end
else
    sol = solveCobraLP(QuickProblem);
end
relSol = sol.full(Order(Presence));
%Obtain fluxes at their boundaries
maxSolved = model.ub(Order(Presence)) == relSol;
minSolved = model.lb(Order(Presence)) == relSol;
if minNorm
    preCompMaxSols(maxSolved) = {sol};
    preCompMinSols(minSolved) = {sol};
end

%Minimise reactions
QuickProblem.osense = 1;
if allowLoops <= 0
    if allowLoops >= -1
        sol = solveCobraMILP(QuickProblem);
    end
else
    sol = solveCobraLP(QuickProblem);
end
relSol = sol.full(Order(Presence));
%Again obtain fluxes at their boundaries
maxSolved = maxSolved | (model.ub(Order(Presence)) == relSol);
minSolved = minSolved | (model.lb(Order(Presence)) == relSol);
%This is only necessary, if we want a min norm.
if minNorm
    preCompMaxSols((model.ub(Order(Presence)) == relSol)) = {sol};
    preCompMinSols((model.lb(Order(Presence)) == relSol)) = {sol};
end
%Restrict the reactions to test only those which are not at their boundariestestFv.
rxnListMin = rxnNameList(~minSolved);
rxnListMax = rxnNameList(~maxSolved);

% generate the loopless problem beforehand instead of during each loop
if allowLoops > 0
    MILPproblem = [];
elseif allowLoops <= 0 && allowLoops >= -1
    MILPproblem = addLoopLawConstraints(LPproblem, model, 1:nRxns, -allowLoops + 1, loopInfo);
else
    % allowLoops <= -2. Use LLCs
    [MILPproblem, loopInfo] = addLoopLawConstraints(LPproblem, model, 1:nRxns, 3, loopInfo);
    rhs0 = MILPproblem.b;
    % no need to regenerate the preprocessing information. They remain unchanged
end

if ~PCT_status || (~exist('parpool') || poolsize == 0)  %aka nothing is active
    
     
    
    if minNorm
        for i = 1:length(rxnNameList)
            
            if allowLoops >= -1
                LPproblem.osense = 1;
                [minFlux(i),Vmin(:,i)] = calcSolForEntry(model,rxnNameList,i,LPproblem,0, ...
                    method, allowLoops,printLevel,minNorm,cpxControl,preCompMinSols{i}, MILPproblem);
                LPproblem.osense = -1;
                [maxFlux(i),Vmax(:,i)] = calcSolForEntry(model,rxnNameList,i,LPproblem,0, ...
                    method, allowLoops,printLevel,minNorm,cpxControl,preCompMaxSols{i}, MILPproblem);
            else
                % apply localized loopless constraints
                i0 = findRxnIDs(model, rxnNameList(i));
                if alwaysLLC || any(loopInfo.rxnInLoops(i0, :))
                    % reset the bounds and rhs in the MILP
                    MILPproblem = restoreOriginalBounds(MILPproblem, rhs0, loopInfo);
                    if ~any(loopInfo.rxnInLoops(i0, :))
                        % apply LLCs only to the always-on set of reactions
                        rxnID = [];
                    else
                        % apply LLCs to the always-on set + objective reaction that is in cycles
                        rxnID = i0;
                    end
                    % update bounds and rhs
                    MILPproblem = updateLLCs(MILPproblem, conCompAlwaysOn, rxnInLoopsAlwaysOn, loopInfo, rxnID, useRxnLink);
                end
                % minimization
                LPproblem.osense = 1;
                if ~alwaysLLC && ~loopInfo.rxnInLoops(i0, 1)
                    % solve as LP is fine if no LLCs are always on and the
                    % reverse direction of the current reaction is not in cycles
                    [minFlux(i),Vmin(:,i)] = calcSolForEntry(model,rxnNameList,i,LPproblem,0, ...
                        method, 1, printLevel,minNorm,cpxControl,preCompMinSols{i}, []);
                else
                    [minFlux(i),Vmin(:,i)] = calcSolForEntry(model,rxnNameList,i,LPproblem,0, ...
                        method, allowLoops, printLevel,minNorm,cpxControl,preCompMinSols{i}, MILPproblem);
                end
                % maximization
                LPproblem.osense = -1;
                if ~alwaysLLC && ~loopInfo.rxnInLoops(i0, 2)
                    % solve as LP is fine if no LLCs are always on and the
                    % forward direction of the current reaction is not in cycles
                    [maxFlux(i),Vmax(:,i)] = calcSolForEntry(model,rxnNameList,i,LPproblem,0, ...
                        method, 1,printLevel,minNorm,cpxControl,preCompMaxSols{i}, []);
                else
                    [maxFlux(i),Vmax(:,i)] = calcSolForEntry(model,rxnNameList,i,LPproblem,0, ...
                        method, allowLoops,printLevel,minNorm,cpxControl,preCompMaxSols{i}, MILPproblem);
                end
            end

            if printLevel == 1 && ~parallelMode
                showprogress(i/length(rxnNameList));
            end
            if printLevel > 1 && ~parallelMode
                fprintf('%4d\t%4.0f\t%10s\t%9.3f\t%9.3f\n',i,100*i/length(rxnNameList),rxnNameList{i},minFlux(i),maxFlux(i));
            end
        end
    else
        %Calc minimums
        mins = -inf*ones(length(rxnListMin),1);
        LPproblem.osense = 1;
        for i = 1:length(rxnListMin)
            if allowLoops >= -1
                [mins(i)] = calcSolForEntry(model,rxnListMin,i,LPproblem,0, method, allowLoops,printLevel,minNorm,cpxControl,[], MILPproblem);
            else
                i0 = findRxnIDs(model, rxnListMin(i));
                if ~alwaysLLC && ~loopInfo.rxnInLoops(i0, 1)
                    % solve as LP is fine if no LLCs are always on and the
                    % reverse direction of the current reaction is not in cycles
                    mins(i) = calcSolForEntry(model,rxnNameList,i,LPproblem,0, ...
                        method, 1, printLevel,minNorm,cpxControl,[], []);
                else
                    % reset the bounds and rhs in the MILP
                    MILPproblem = restoreOriginalBounds(MILPproblem, rhs0, loopInfo);
                    rxnID = [];
                    if loopInfo.rxnInLoops(i0, 1)  % if the reverse direction of rxn i0 is in cycles
                        % apply LLCs to the objective reaction
                        rxnID = i0;
                    end
                    % update bounds and rhs
                    MILPproblem = updateLLCs(MILPproblem, conCompAlwaysOn, rxnInLoopsAlwaysOn, loopInfo, rxnID, useRxnLink);
                    mins(i) = calcSolForEntry(model,rxnNameList,i,LPproblem,0, ...
                        method, allowLoops, printLevel,minNorm,cpxControl,[], MILPproblem);
                end
            end
        end
        [minFluxPres,minFluxOrder] = ismember(rxnListMin,rxnNameList);
        minFlux(minFluxOrder(minFluxPres)) = mins;   
        %calc maximiums
        maxs = inf*ones(length(rxnListMax),1);
        LPproblem.osense = -1;
        for i = 1:length(rxnListMax)
            if allowLoops >= -1
                [maxs(i)] = calcSolForEntry(model,rxnListMax,i,LPproblem,0, method, allowLoops,printLevel,minNorm,cpxControl,[], MILPproblem);
            else
                i0 = findRxnIDs(model, rxnListMax(i));
                if ~alwaysLLC && ~loopInfo.rxnInLoops(i0, 2)
                    % solve as LP is fine if no LLCs are always on and the
                    % forward direction of the current reaction is not in cycles
                    maxs(i) = calcSolForEntry(model,rxnNameList,i,LPproblem,0, ...
                        method, 1, printLevel,minNorm,cpxControl,[], []);
                else
                    % reset the bounds and rhs in the MILP
                    MILPproblem = restoreOriginalBounds(MILPproblem, rhs0, loopInfo);
                    rxnID = [];
                    if loopInfo.rxnInLoops(i0, 2)  % if the forward direction of rxn i0 is in cycles
                        % apply LLCs to the objective reaction
                        rxnID = i0;
                    end
                    % update bounds and rhs
                    MILPproblem = updateLLCs(MILPproblem, conCompAlwaysOn, rxnInLoopsAlwaysOn, loopInfo, rxnID, useRxnLink);
                    maxs(i) = calcSolForEntry(model,rxnNameList,i,LPproblem,0, ...
                        method, allowLoops, printLevel,minNorm,cpxControl,[], MILPproblem);
                end
            end
        end
        [maxFluxPres,maxFluxOrder] = ismember(rxnListMax,rxnNameList);
        maxFlux(maxFluxOrder(maxFluxPres)) = maxs; 
    end
else % parallel job.  pretty much does the same thing.

    global CBT_LP_SOLVER;
    global CBT_MILP_SOLVER;
    global CBT_QP_SOLVER;
    lpsolver = CBT_LP_SOLVER;
    qpsolver = CBT_QP_SOLVER;
    milpsolver = CBT_MILP_SOLVER;    
    if minNorm        
        parfor i = 1:length(rxnNameList)
            changeCobraSolver(qpsolver,'QP',0,1);
            changeCobraSolver(lpsolver,'LP',0,1);
            changeCobraSolver(milpsolver,'MILP',0,1);
            
            parLPproblem = LPproblem;
            parMILPproblem = MILPproblem;
            if allowLoops >= -1    
                parLPproblem.osense = 1;
                [minFlux(i),Vmin(:,i)] = calcSolForEntry(model,rxnNameList,i,parLPproblem,1, method, allowLoops,printLevel,minNorm,cpxControl,preCompMinSols{i}, parMILPproblem);
                parLPproblem.osense = -1;
                [maxFlux(i),Vmax(:,i)] = calcSolForEntry(model,rxnNameList,i,parLPproblem,1, method, allowLoops,printLevel,minNorm,cpxControl,preCompMaxSols{i}, parMILPproblem);
            else
                % apply localized loopless constraints
                i0 = findRxnIDs(model, rxnNameList(i));
                if alwaysLLC || any(loopInfo.rxnInLoops(i0, :))
                    if ~any(loopInfo.rxnInLoops(i0, :))
                        % apply LLCs only to the always-on set of reactions
                        rxnID = [];
                    else
                        % apply LLCs to the always-on set + objective reaction that is in cycles
                        rxnID = i0;
                    end
                    % update bounds and rhs
                    parMILPproblem = updateLLCs(parMILPproblem, conCompAlwaysOn, rxnInLoopsAlwaysOn, loopInfo, rxnID, useRxnLink);
                end
                % minimization
                parLPproblem.osense = 1;
                if ~alwaysLLC && ~loopInfo.rxnInLoops(i0, 1)
                    % solve as LP is fine if no LLCs are always on and the
                    % reverse direction of the current reaction is not in cycles
                    [minFlux(i),Vmin(:,i)] = calcSolForEntry(model,rxnNameList,i, parLPproblem,0, ...
                        method, 1, printLevel,minNorm,cpxControl,preCompMinSols{i}, []);
                else
                    [minFlux(i),Vmin(:,i)] = calcSolForEntry(model,rxnNameList,i, parLPproblem,0, ...
                        method, allowLoops, printLevel,minNorm,cpxControl,preCompMinSols{i}, parMILPproblem);
                end
                % maximization
                parLPproblem.osense = -1;
                if ~alwaysLLC && ~loopInfo.rxnInLoops(i0, 2)
                    % solve as LP is fine if no LLCs are always on and the
                    % forward direction of the current reaction is not in cycles
                    [maxFlux(i),Vmax(:,i)] = calcSolForEntry(model,rxnNameList,i, parLPproblem,0, ...
                        method, 1,printLevel,minNorm,cpxControl,preCompMaxSols{i}, []);
                else
                    [maxFlux(i),Vmax(:,i)] = calcSolForEntry(model,rxnNameList,i, parLPproblem,0, ...
                        method, allowLoops,printLevel,minNorm,cpxControl,preCompMaxSols{i}, parMILPproblem);
                end
            end
        end
    else
        mins = -inf*ones(length(rxnListMin),1);
        LPproblem.osense = 1;
        parfor i = 1:length(rxnListMin)    
            changeCobraSolver(qpsolver,'QP',0,1);
            changeCobraSolver(lpsolver,'LP',0,1);
            parLPproblem = LPproblem;
            parMILPproblem = MILPproblem;
            if allowLoops >= -1
                [mins(i)] = calcSolForEntry(model,rxnListMin,i,parLPproblem,1, method, allowLoops,printLevel,minNorm,cpxControl,[], parMILPproblem);
            else
                i0 = findRxnIDs(model, rxnListMin(i));
                if ~alwaysLLC && ~loopInfo.rxnInLoops(i0, 1)
                    % solve as LP is fine if no LLCs are always on and the
                    % reverse direction of the current reaction is not in cycles
                    mins(i) = calcSolForEntry(model,rxnNameList,i,parLPproblem,0, ...
                        method, 1, printLevel,minNorm,cpxControl,[], []);
                else
                    if loopInfo.rxnInLoops(i0, 1)  % if the reverse direction of rxn i0 is in cycles
                        % apply LLCs only to the always-on set of reactions
                        rxnID = [];
                    else
                        % apply LLCs to the always-on set + objective reaction that is in cycles
                        rxnID = i0;
                    end
                    % update bounds and rhs
                    parMILPproblem = updateLLCs(parMILPproblem, conCompAlwaysOn, rxnInLoopsAlwaysOn, loopInfo, rxnID, useRxnLink);
                    mins(i) = calcSolForEntry(model,rxnNameList,i,parLPproblem,0, ...
                        method, allowLoops, printLevel,minNorm,cpxControl,[], parMILPproblem);
                end
            end
        end
        [minFluxPres,minFluxOrder] = ismember(rxnListMin,rxnNameList);
        minFlux(minFluxOrder(minFluxPres)) = mins;   
        %calc maximiums
        maxs = inf*ones(length(rxnListMax),1);
        LPproblem.osense = -1;
        parfor i = 1:length(rxnListMax)        
            changeCobraSolver(qpsolver,'QP',0,1);
            changeCobraSolver(lpsolver,'LP',0,1);
            parLPproblem = LPproblem;
            parMILPproblem = MILPproblem;
            if allowLoops >= -1
                [maxs(i)] = calcSolForEntry(model,rxnListMax,i,parLPproblem,1, method, allowLoops,printLevel,minNorm,cpxControl,[], parMILPproblem);
            else
                i0 = findRxnIDs(model, rxnListMax(i));
                if ~alwaysLLC && ~loopInfo.rxnInLoops(i0, 2)
                    % solve as LP is fine if no LLCs are always on and the
                    % forward direction of the current reaction is not in cycles
                    maxs(i) = calcSolForEntry(model,rxnNameList,i, parLPproblem,0, ...
                        method, 1, printLevel,minNorm,cpxControl,[], []);
                else
                    if loopInfo.rxnInLoops(i0, 2)  % if the forward direction of rxn i0 is in cycles
                        % apply LLCs only to the always-on set of reactions
                        rxnID = [];
                    else
                        % apply LLCs to the always-on set + objective reaction that is in cycles
                        rxnID = i0;
                    end
                    % update bounds and rhs
                    parMILPproblem = updateLLCs(parMILPproblem, conCompAlwaysOn, rxnInLoopsAlwaysOn, loopInfo, rxnID, useRxnLink);
                    maxs(i) = calcSolForEntry(model,rxnNameList,i, parLPproblem,0, ...
                        method, allowLoops, printLevel,minNorm,cpxControl,[], parMILPproblem);
                end
            end
        end
        [maxFluxPres,maxFluxOrder] = ismember(rxnListMax,rxnNameList);
        maxFlux(maxFluxOrder(maxFluxPres)) = maxs;         
    end
end

maxFlux = columnVector(maxFlux);
minFlux = columnVector(minFlux);

    
end

function [Flux,V] = calcSolForEntry(model,rxnNameList,i,LPproblem,parallelMode, method, allowLoops, printLevel, minNorm, cpxControl, sol, MILPproblem)

    %get Number of reactions
    nRxns = numel(model.rxns);
    %Set the correct objective
    LPproblem.c = double(ismember(model.rxns,rxnNameList{i}));
    if isempty(sol)
        if printLevel == 1 && ~parallelMode
            fprintf('iteration %d.\n', i);
        end
        % do LP always
        if allowLoops > 0
            LPsolution = solveCobraLP(LPproblem, cpxControl);
        else
            MILPproblem.osense = LPproblem.osense;
            MILPproblem.c(1:nRxns) = LPproblem.c;
            LPsolution = solveCobraMILP(MILPproblem);
        end
        % take the maximum flux from the flux vector, not from the obj -Ronan
        % A solution is possible, so the only problem should be if its
        % unbounded and if it is unbounded, the max flux is infinity.
        if LPsolution.stat == 2
            Flux = -LPProblem.osense * inf;
        else
            Flux = getObjectiveFlux(LPsolution, LPproblem);
        end
    else
        LPsolution = sol;
        Flux = getObjectiveFlux(LPsolution, LPproblem);        
    end
    % minimise the Euclidean norm of the optimal flux vector to remove loops -Ronan
    if minNorm == 1
        V = getMinNorm(LPproblem, LPsolution, nRxns, Flux, model, method);
    end
end


function V = getMinNorm(LPproblem,LPsolution,nRxns,cFlux, model, method)
% get the Flux distribution for the specified min norm.

    if strcmp(method, '2-norm')
        QPproblem=LPproblem;
        QPproblem.lb(LPproblem.c~=0) = cFlux - 1e-12;
        QPproblem.ub(LPproblem.c~=0) = cFlux + 1e12;
        QPproblem.c(:)=0;
        %Minimise Euclidean norm using quadratic programming
        QPproblem.F = speye(nRxns,nRxns);
        QPproblem.osense = 1;
        %quadratic optimization
        solution = solveCobraQP(QPproblem);
        V=solution.full(1:nRxns,1);
    elseif strcmp(method, '1-norm')
        vSparse = sparseFBA(LPproblem, 'min', 0, 0, 'l1');
        V = vSparse;
    elseif strcmp(method, '0-norm')
        vSparse = sparseFBA(LPproblem, 'min', 0, 0);
        V = vSparse;
    elseif strcmp(method, 'FBA')
        V=LPsolution.full(1:nRxns);
    elseif strcmp(method, 'minOrigSol')
        LPproblemMOMA = LPproblem;
        LPproblemMOMA=rmfield(LPproblemMOMA, 'csense');
        LPproblemMOMA.A = model.S;
        LPproblemMOMA.S = LPproblemMOMA.A;
        LPproblemMOMA.b = model.b;
        LPproblemMOMA.lb(LPproblem.c~=0) = cFlux - 1e-12;
        LPproblemMOMA.ub(LPproblem.c~=0) = cFlux + 1e-12;
        LPproblemMOMA.rxns = model.rxns;
        momaSolution = linearMOMA(model,LPproblemMOMA);
        V=momaSolution.x;
    end
end


function flux = getObjectiveFlux(LPsolution,LPproblem)
% Determine the current flux based on an LPsolution, the original LPproblem
% The LPproblem is used to retrieve the current objective position.
% min indicates, whether the minimum or maximum is requested, the
% upper/lower bounds are used, if the value is exceeding them

    Index = LPproblem.c~=0;
    if LPsolution.full(Index)<LPproblem.lb(Index) %takes out tolerance issues
        flux = LPproblem.lb(Index);
    elseif LPsolution.full(Index)>LPproblem.ub(Index)
        flux = LPproblem.ub(Index);
    else
        flux = LPsolution.full(Index);
    end
end

function [alwaysLLC, rxnInLoopsAlwaysOn, conCompAlwaysOn, x0] = preprocessLLC(LPproblem, ...
    model, nRxns, rxnInLoops, osenseStr, conComp, printLevel)
if strcmp(osenseStr, 'min')
    model.c = -model.c;
end
% determine the set of reactions for which LLCs are always required
% condition I in Prop. 2 in Chan et al., 2017
cond1 = rxnInLoops(:, 2) & model.c > 0;
% condition II in the paper in Prop. 2 in Chan et al., 2017 
cond2 = rxnInLoops(:, 1) & model.c < 0;
% condition III in the paper in Prop. 2 in Chan et al., 2017
[cond3A1, cond3A2, cond3B] = deal(false(nRxns, 1));
for i = (size(model.S, 1) + 1):size(LPproblem.A, 1)
    % for constraint p with sum(a_pj * v_j) <= b_p
    if ~strcmp(LPproblem.csense(i), 'G')  % '<=' or '=' constraint
        % if reaction j has its forward direction in cycles and a_pj < 0
        cond3A1 = cond3A1 | (rxnInLoops(:, 2) & LPproblem.A(i, 1:nRxns)' < 0);
        % if reaction j has its reverse direction in cycles and a_pj > 0
        cond3A2 = cond3A2 | (rxnInLoops(:, 1) & LPproblem.A(i, 1:nRxns)' > 0);
        % if the constraint involves 2 or more reactinos or RHS < 0
        cond3B = cond3B | (nnz(LPproblem.A(i, 1:nRxns)) > 1 | LPproblem.b(i) < 0);
    end
    if ~strcmp(LPproblem.csense(i), 'L')  % '>=' or '=' constraint
        cond3A1 = cond3A1 | (rxnInLoops(:, 2) & LPproblem.A(i, 1:nRxns)' > 0);
        cond3A2 = con3A2 | (rxnInLoops(:, 1) & LPproblem.A(i, 1:nRxns)' < 0);
        cond3B = cond3B | (nnz(LPproblem.A(i, 1:nRxns)) > 1 | LPproblem.b(i) > 0);
    end
end
% reactions satisfying (3A1 or 3A2) and 3B
cond3 = (cond3A1 | cond3A2) & cond3B;
% condition III for bound constraints can be simplified as follows:
cond3 = cond3 | (model.lb > 0 & rxnInLoops(:, 2)) | (model.ub < 0 & rxnInLoops(:, 1));
% reactions that are required to be constrained by loopless constraints all the time
rxnInLoopsAlwaysOn = cond1 | cond2 | cond3;
% LLCs are always required if the set is non-empty
alwaysLLC = any(rxnInLoopsAlwaysOn);
% the corresponding set of reactions in the same connected components as
% the always-on reactions
conCompAlwaysOn = false(max(conComp), 1);
conCompAlwaysOn(conComp(rxnInLoopsAlwaysOn)) = true;
if printLevel 
    fprintf('Reactions in internal nullspace can be divided into %d connected components.\n', max(conComp))
end

% get an initial feasible and loopless solution in case MipStart is needed
model2 = model;
model2.lb = model2.lb(1:size(model2.S, 2));
model2.ub = model2.ub(1:size(model2.S, 2));
model2.c = zeros(size(model2.S, 2), 1);
model2.b = zeros(size(model2.S, 1), 1);
sFeas = optimizeCbModel(model2, 'max', 'one');
x0 = sFeas.x;
clear sFeas model2
end

function MILPproblemLLC = updateLLCs(MILPproblemLLC, conCompAlwaysOn, rxnInLoopsAlwaysOn, llcInfo, rxnID, useRxnLink)
% apply LLCs by relaxing constraints and pre-assign values to variables
if nargin < 2
    rxnID = [];
end
conCompOn = conCompAlwaysOn;
conCompOn(llcInfo.conComp(rxnID)) = true;

bigM = inf;
if ~useRxnLink
    % use connections from nullspace
    for jCon = 1:numel(conCompOn)
        if ~conCompOn(jCon)
            % relax constraints not affecting optimality and feasibility
            MILPproblemLLC.b(llcInfo.con.vU(llcInfo.rxnInLoopIds(llcInfo.conComp == jCon))) = bigM;
            MILPproblemLLC.b(llcInfo.con.gU(llcInfo.rxnInLoopIds(llcInfo.conComp == jCon))) = bigM;
            MILPproblemLLC.b(llcInfo.con.vL(llcInfo.rxnInLoopIds(llcInfo.conComp == jCon))) = -bigM;
            MILPproblemLLC.b(llcInfo.con.gL(llcInfo.rxnInLoopIds(llcInfo.conComp == jCon))) = bigM;
            % fix variables not affecting optimality and feasibility
            MILPproblemLLC.lb(llcInfo.var.g(llcInfo.rxnInLoopIds(llcInfo.conComp == jCon))) = 0;
            MILPproblemLLC.ub(llcInfo.var.g(llcInfo.rxnInLoopIds(llcInfo.conComp == jCon))) = 0;
            MILPproblemLLC.ub(llcInfo.var.z(llcInfo.rxnInLoopIds(llcInfo.conComp == jCon))) = 0;
        end
    end
else
    % use connections from EFMs
    rxnOn = rxnInLoopsAlwaysOn;
    rxnOn(rxnID) = true;
    % reactions in cycles not sharing EFMs with the current rxn and
    % not one of the reactions required to have no flux through cycles
    id = ~any(llcInfo.rxnLink(rxnOn, :), 1)' & any(llcInfo.rxnInLoops, 2);
%     for iRxn = find(rxnOn(:))'
        
%         id = llcInfo.rxnLink(iRxn, :)' ~= 0
%     end
    MILPproblemLLC.b(llcInfo.con.vU(llcInfo.rxnInLoopIds(id))) = bigM;
    MILPproblemLLC.b(llcInfo.con.gU(llcInfo.rxnInLoopIds(id))) = bigM;
    MILPproblemLLC.b(llcInfo.con.vL(llcInfo.rxnInLoopIds(id))) = -bigM;
    MILPproblemLLC.b(llcInfo.con.gL(llcInfo.rxnInLoopIds(id))) = bigM;

    % pre-determine variables not connected to the reaction for FVA
    % except reactions required to be always constrained
    rxnKeep = llcInfo.conComp == 0;
    for jCon = 1:numel(conCompOn)
        if conCompOn(jCon)
            rxnKeep(llcInfo.conComp == jCon) = true;
        end
    end
    MILPproblemLLC.lb(llcInfo.var.g(llcInfo.rxnInLoopIds(~rxnKeep))) = 0;
    MILPproblemLLC.ub(llcInfo.var.g(llcInfo.rxnInLoopIds(~rxnKeep))) = 0;
    MILPproblemLLC.ub(llcInfo.var.z(llcInfo.rxnInLoopIds(~rxnKeep))) = 0;
end
end

function MILPproblemLLC = restoreOriginalBounds(MILPproblemLLC, rhs0, llcInfo)
    MILPproblemLLC.b = rhs0;
    MILPproblemLLC.ub(llcInfo.var.z) = 1;
    MILPproblemLLC.ub(llcInfo.var.g) = llcInfo.BDg;
    MILPproblemLLC.lb(llcInfo.var.g) = -llcInfo.BDg;
end
