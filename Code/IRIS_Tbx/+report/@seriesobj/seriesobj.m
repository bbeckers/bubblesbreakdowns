classdef seriesobj < report.genericobj & report.condformatobj
    
    properties
        data = {tseries()};
    end
    
    methods
        
        function This = seriesobj(varargin)
            This = This@report.genericobj(varargin{:});
            This = This@report.condformatobj();
            This.childof = {'graph','table'};
            This.default = [This.default,{ ...
                'autodata',{},@(x) isempty(x) ...
                || isfunc(x) ...
                || iscell(x),true, ...
                'colstruct',[],@(x) isempty(x) ...
                || report.genericobj.validatecolstruct(x),true, ...
                'condformat',[], ...
                @(x) isempty(x) || ( ...
                isstruct(x) ...
                && isfield(x,'test') && isfield(x,'format') ...
                && iscellstr({x.test}) && iscellstr({x.format}) ), ...
                true, ...
                'decimal',NaN,@isnumericscalar,true, ...
                'format','%.2f',@ischar,true, ...
                'highlight',[],@(x) isnumeric(x) || isfunc(x),true, ...
                'legend',Inf,@(x) isempty(x) ...
                || (isnumericscalar(x) && (isnan(x) || isinf(x))) ...
                || iscellstr(x) || ischar(x),false, ...
                ... 'lhs',Inf,@(x) isnumericscalar(x) && ...
                ... (isequal(x,Inf) || (x == round(x) && x >= 1)), ...
                ... true, ...
                'marks',{},@(x) isempty(x) || iscellstr(x),true, ...
                'inf','\ensuremath{\infty}',@ischar,true, ...
                'nan','\ensuremath{\cdots}',@ischar,true, ...
                'plotfunc',@plot,@(x) isfunc(x) ...
                && any(strcmpi(char(x),{'plot','area','bar','conbar', ...
                'barcon','stem','plotcmp','plotpred'})), ...
                true, ...
                'plotoptions',{}, ...
                @(x) iscell(x) && iscellstr(x(1:2:end)), ...
                true, ...
                'purezero','',@ischar,true, ...
                'printedzero','',@ischar,true, ...
                'separator','',@ischar,false, ...
                'showmarks',true,@islogical,true, ...
                'units','',@ischar,true, ...
                'yaxis,yaxislocation','left',@(x) any(strcmpi(x,{'left','right'})),false, ...
                }];
        end
        
        function [This,varargin] = specargin(This,varargin)
            if ~isempty(varargin)
                if istseries(varargin{1}) || iscell(varargin{1})
                    This.data = varargin{1};
                    if istseries(This.data)
                        This.data = {This.data};
                    end
                end
                varargin(1) = [];
            end
        end
        
        function This = setoptions(This,varargin)
            This = setoptions@report.genericobj(This,varargin{:});
            This.options.marks = This.options.marks(:).';
            if ~isempty(This.options.autodata)
                if isa(This.options.autodata,'function_handle')
                    This.options.autodata = {This.options.autodata};
                end
                This = autodata(This);
            end
            if ~isnan(This.options.decimal)
                This.options.format = sprintf('%%.%gf', ...
                    round(abs(This.options.decimal)));
            end
            This = assign(This,This.options.condformat);
            if ~This.options.showmarks
                This.options.marks = {};
            end
            if strcmp(char(This.options.plotfunc),'conbar')
                This.options.plotfunc = @barcon;
            end
        end
        
        varargout = plot(varargin)
        varargout = latexonerow(varargin)
        varargout = latexdata(varargin)
        varargout = autodata(varargin)
        varargout = getdata(varargin)
        
    end
    
    methods (Access=protected,Hidden)
        varargout = mylegend(varargin)
        varargout = speclatexcode(varargin)
    end

end