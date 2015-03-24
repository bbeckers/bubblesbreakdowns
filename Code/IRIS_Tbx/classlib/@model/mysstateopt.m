function Opt = mysstateopt(This,Mode,varargin)
% mysstateopt  [Not a public function] Prepare steady-state solver options.
%
% Backend IRIS function.
% No help provided.

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

if length(varargin) == 1 && isequal(varargin{1},false)
    Opt = false;
    return
end

if length(varargin) == 1 && isequal(varargin{1},true)
    varargin(1) = [];
end

% `Mode` is either `'verbose'` (direct calls to `model/sstate`) or
% `'silent'`; the mode determines the default values for `'display='` and
% `'warning='`.
Opt = passvalopt(['model.mysstate',Mode],varargin{:});

%--------------------------------------------------------------------------

if This.linear
    % Linear sstate solver
    %----------------------
    % No need to pre-process any options for the linear sstate solver.
else
    % Non-linear sstate solver
    %--------------------------
    Opt = xxBlocks(This,Opt);
    Opt = xxDisplayOpt(This,Opt);
    Opt = xxOptimOpt(This,Opt);
end

end

% Subfunctions.

%**************************************************************************
function Opt = xxDisplayOpt(This,Opt) %#ok<INUSL>
if islogical(Opt.display)
    if Opt.display
        Opt.display = 'iter';
    else
        Opt.display = 'off';
    end
end
end % xxDisplayOpt().

%**************************************************************************
function Opt = xxOptimOpt(This,Opt) %#ok<INUSL>
% Use Levenberg-Marquardt because it can handle underdetermined systems.
oo = Opt.optimset;
if ~isempty(oo)
    oo(1:2:end) = regexprep(oo(1:2:end),'[^\w]','');
end
Opt.optimset = optimset( ...
    'display',Opt.display, ...
    'maxiter',Opt.maxiter, ...
    'maxfunevals',Opt.maxfunevals,  ...
    'tolx',Opt.tolx, ...
    'tolfun',Opt.tolfun, ...
    'algorithm','levenberg-marquardt', ...
    oo{:});
end % xxOptimOpt().

%**************************************************************************
function Opt = xxBlocks(This,Opt)

% Process fix options first.
fixL = [];
fixG = [];
doFixOpt();

% Swap nametype of exogenised variables and endogenised parameters.
isSwap = ~isempty(Opt.endogenise) || ~isempty(Opt.exogenise);
if isSwap
    This = mysstateswap(This,Opt);
end

% Run BLAZER if it has not been run yet or if user requested
% exogenise/endogenise.
if Opt.blocks && (isempty(This.nameblk) || isempty(This.eqtnblk) || isSwap)
    This = myblazer(This);
end

% Prepare blocks of equations/names.
if Opt.blocks
    nameBlkL = This.nameblk;
    nameBlkG = This.nameblk;
    eqtnBlk = This.eqtnblk;
else
    % If `'blocks=' false`, prepare two blocks:
    % # transition equations;
    % # measurement equations.
    nameBlkL = cell(1,2);
    nameBlkG = cell(1,2);
    % All transition equations and variables.
    eqtnBlk = cell(1,2);
    nameBlkL{1} = find(This.nametype == 2);
    nameBlkG{1} = find(This.nametype == 2);
    eqtnBlk{1} = find(This.eqtntype == 2);
    % All measurement equations and variables.
    nameBlkL{2} = find(This.nametype == 1);
    nameBlkG{2} = find(This.nametype == 1);
    eqtnBlk{2} = find(This.eqtntype == 1);
end

nBlk = length(nameBlkL);
blkFunc = cell(1,nBlk);
% Remove variables fixed by the user.
% Prepare function handles to evaluate individual equation blocks.
for ii = 1 : nBlk
    % Exclude fixed levels and growth rates from the list of optimised
    % names.
    nameBlkL{ii} = setdiff(nameBlkL{ii},fixL);
    nameBlkL{ii} = setdiff(nameBlkL{ii},This.Refresh);
    nameBlkG{ii} = setdiff(nameBlkG{ii},fixG);
    nameBlkG{ii} = setdiff(nameBlkG{ii},This.Refresh);
    if isempty(nameBlkL{ii}) && isempty(nameBlkG{ii})
        continue
    end
    % Create an anonymous function handle for each block.
    eqtn = This.eqtnS(eqtnBlk{ii});
    eqtn = strrep(eqtn,'exp?','exp');
    % Replace log(exp(x(...))) with x(...). This helps a lot.
    eqtn = regexprep(eqtn,'log\(exp\(x\((\d+)\)\)\)','x($1)');
    % Create a function handle used to evaluate each block of
    % equations.
    blkFunc{ii} = eval(['@(x,dx) [',eqtn{:},']']);
end

% Index of level and growth variables endogenous in sstate calculation.
endogLInx = false(size(This.name));
endogLInx([nameBlkL{:}]) = true;
endogGInx = false(size(This.name));
endogGInx([nameBlkG{:}]) = true;

% Index of level variables that will be always set to zero.
zeroLInx = false(size(This.name));
zeroLInx(This.nametype == 3) = true;
zeroGInx = false(size(This.name));
zeroGInx(This.nametype == 3) = true;

Opt.fixL = fixL;
Opt.fixG = fixG;
Opt.nameBlkL = nameBlkL;
Opt.nameBlkG = nameBlkG;
Opt.eqtnBlk = eqtnBlk;
Opt.blkFunc = blkFunc;
Opt.endogLInx = endogLInx;
Opt.endogGInx = endogGInx;
Opt.zeroLInx = zeroLInx;
Opt.zeroGInx = zeroGInx;

    function doFixOpt()
        % Process the fix, fixallbut, fixlevel, fixlevelallbut, fixgrowth,
        % and fixgrowthallbut options. All the user-supply information is
        % combined into fixlevel and fixgrowth.
        canBeFixed = This.nametype <= 2 | This.nametype == 4;
        list = {'fix','fixlevel','fixgrowth'};
        for i = 1 : length(list)
            fix = list{i};
            fixAllBut = [fix,'allbut'];
            
            % Convert charlist to cellstr.
            if ischar(Opt.(fix)) ...
                    && ~isempty(Opt.(fix))
                Opt.(fix) = regexp(Opt.(fix),'\w+','match');
            end
            if ischar(Opt.(fixAllBut)) ...
                    && ~isempty(Opt.(fixAllBut))
                Opt.(fixAllBut) = ...
                    regexp(Opt.(fixAllBut),'\w+','match');
            end
            
            % Convert fixAllBut to fix.
            if ~isempty(Opt.(fixAllBut))
                Opt.(fix) = ...
                    setdiff(This.name(canBeFixed),Opt.(fixAllBut));
            end
            
            if ~isempty(Opt.(fix))
                fixpos = mynameposition(This,Opt.(fix));
                validate = ~isnan(fixpos);
                if all(validate)
                    validate = This.nametype(fixpos) <= 2 ...
                        | This.nametype(fixpos) == 4;
                end
                if any(~validate)
                    utils.error('model', ...
                        'Cannot fix this name: ''%s''.', ...
                        Opt.(fix){~validate});
                end
                Opt.(fix) = fixpos;
            else
                Opt.(fix) = [];
            end
        end
        
        % Add the positions of optimal policy multipliers to the list of fixed
        % variables. The level and growth of multipliers will be set to zero in the
        % main loop.
        if Opt.zeromultipliers
            Opt.fix = union(Opt.fix,find(This.multiplier));
        end
        
        fixL = union(Opt.fix,Opt.fixlevel);
        if Opt.growth
            fixG = union(Opt.fix,Opt.fixgrowth);
        else
            fixG = find(canBeFixed);
        end
    end % doOpt2Fix().

end % xxBlocks().