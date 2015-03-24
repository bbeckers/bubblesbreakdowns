function [Code,Labels,Export,Subs,Comment] ...
    = readcode(FileList,Params,Labels,Export,ParentFile,Opt)
% readcode  [Not a public function] Preparser master file.
%
% Backend IRIS function.
% No help provided.

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

if isempty(Params)
    Params = struct([]);
end

if isempty(Labels)
    Labels = {};
end

if isempty(Export)
    Export = struct('filename',{},'content',{});
end

if isempty(ParentFile)
    ParentFile = '';
end

if isempty(Opt.removecomments)
    Opt.removecomments = {};
end

%--------------------------------------------------------------------------

Code = '';
Comment = '';
Subs = struct();

if ischar(FileList)
    FileList = {FileList};
end

fileStr = '';
nFileList = length(FileList);
for i = 1 : nFileList
    Code = [Code,sprintf('\n'),file2char(FileList{i})]; %#ok<AGROW>
    fileStr = [fileStr,FileList{i}]; %#ok<AGROW>
    if i < nFileList
        fileStr = [fileStr,' & ']; %#ok<AGROW>
    end
end

fileStr = sprintf('<a href="matlab:edit %s">%s</a>',fileStr,fileStr);
if ~isempty(ParentFile)
    fileStr = [ParentFile,' > ',fileStr];
end
errorParsing = ['Error preparsing file(s) ',fileStr,'. '];

% Convert end of lines.
Code = strfun.converteols(Code);

% Check if there is an initial %% comment line that will be used as comment
% in model objects.
tokens = regexp(Code,'^\s*%%([^\n])+','tokens','once');
if ~isempty(tokens)
    Comment = strtrim(tokens{1});
end

% Read quoted strings 'xxx' and "xxx" and replace them with #(00).
% The quoted strings must be contained at one line.
[Code,Labels] = preparser.protectlabels(Code,Labels);

% Remove standard line and block comments.
Code = strfun.removecomments(Code,Opt.removecomments{:});

% Remove triple exclamation point !!!.
% This mark is meant to be used to highlight some bits of the code.
Code = strrep(Code,'!!!','');

% Characters beyond char(highcharcode) not allowed except comments.
% Default is 1999.
charCap = irisget('highcharcode');
if any(Code > char(charCap))
    utils.error('preparser',[errorParsing, ...
        'The file contains characters beyond char(%g).'],charCap);
end

% Replace @keywords with !keywords. This is for backward compatibility
% only.
Code = xxPrefix(Code);

% Discard everything after !stop.
pos = strfind(Code,'!stop');
if ~isempty(pos)
    Code = Code(1:pos(1)-1);
end

% Add a white space at the end.
Code = [Code,sprintf('\n')];

% Execute control commands
%--------------------------
% Evaluate and expand the following control commands:
% * !if .. !else .. !elseif .. !end
% * !for .. !do .. !end
% * !switch .. !case ... !otherwise .. !end
% * !export .. !end
[Code,Labels,Export] ...
    = preparser.controls(Code,Params,errorParsing,Labels,Export);

% Import external files
%-----------------------

[Code,Labels,Export] = xxImport(Code, ...
    Params,Labels,Export,fileStr,Opt);

% Expand pseudofunctions
%------------------------

[Code,invalid] = preparser.pseudofunc(Code);
if ~isempty(invalid)
    invalid = xxFormatError(invalid,Labels);
    utils.error('preparser',[errorParsing, ...
        'Invalid pseudofunction: ''%s''.'], ...
        invalid{:});
end

% Expand substitutions
%----------------------

% Expand substitutions in the top file after all imports have been done.
if isempty(ParentFile)
    [Code,Subs,leftover,multiple,undef] ...
        = preparser.substitute(Code);
    if ~isempty(leftover)
        leftover = xxFormatError(leftover,Labels);
        utils.error('preparser',[errorParsing, ...
            'There is a leftover code in a substitutions block: ''%s''.'], ...
            leftover{:});
    end
    if ~isempty(multiple)
        multiple = xxFormatError(multiple,Labels);
        utils.error('preparser',[errorParsing, ...
            'This substitution name is defined more than once: ''%s''.'], ...
            multiple{:});
    end
    if ~isempty(undef)
        utils.error('preparser',[errorParsing, ...
            'This is an unresolved or undefined substitution: ''%s''.'], ...
            undef{:});
    end
end

% Remove leading and trailing empty lines.
Code = strfun.removeltel(Code);

end

% Subfunctions.

%**************************************************************************
function Code = xxPrefix(Code)
% doprefix  Replace @keywords with !keywords.

Code = regexprep(Code,'@([a-z]\w*)','!$1');

end % xxPrefix().

%**************************************************************************
function [Code,Labels,Export] ...
    = xxImport(Code,Params,Labels,Export,ParentFile,Opt)

% doimport  Import external file.
% Call import/include/input files with replacement.
pattern = '(?:!include|!input|!import)\((.*?)\)';
while true
    [tokens,start,finish] = ...
        regexp(Code,pattern,'tokens','start','end','once');
    if isempty(tokens)
        break
    end
    fname = strtrim(tokens{1});
    if ~isempty(fname)
        [impcode,Labels,Export] = preparser.readcode(fname, ...
            Params,Labels,Export,ParentFile,Opt);
        Code = [Code(1:start-1),impcode,Code(finish+1:end)];
    end
end

end % xxImport().

%**************************************************************************
function C = xxFormatError(C,Labels)

if iscell(C)
    for i = 1 : length(C)
        C{i} = xxFormatError(C{i},Labels);
    end
    return
end

%C = xxLabelsBack(C,Labels);
C = preparser.labelsback(C,Labels);
C = regexprep(C,'\s+',' ');
C = strtrim(C);
C = strfun.maxdisp(C,40);

end % xxFormatError().