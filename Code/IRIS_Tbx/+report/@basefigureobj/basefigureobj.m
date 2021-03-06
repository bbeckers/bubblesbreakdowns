classdef basefigureobj < report.tabularobj
    
    properties
        handle = [];
    end
    
    methods
        
        function This = basefigureobj(varargin)
            This = This@report.tabularobj(varargin{:});
            This.childof = {'report','align'};
            This.default = [This.default,{ ...
                'close',true,@islogicalscalar,true, ...
                'figureoptions',{}, ...
                @(x) iscell(x) && iscellstr(x(1:2:end)), ...
                true, ...
                'figurescale','auto', ...
                @(x) isnumericscalar(x) || strcmpi(x,'auto'), ...
                true, ...
                'papertype','usletter', ...
                @(x) any(strcmpi(x,{'usletter','uslegal','A4'})), ...
                true, ...
                'subplot','auto', ...
                @(x) (isnumeric(x) && numel(x) == 2) || strcmpi(x,'auto'), ...
                true, ...
                'separator','\medskip\par',@ischar,true, ...
                'style',[],@(x) isempty(x) || isstruct(x),true, ...
                'typeface','',@ischar,false, ...                
                'visible',false,@islogical,true, ...
                }];
        end
        
        % Process class-specific input arguments.
        function [This,varargin] = specargin(This,varargin)
        end
        
        function This = setoptions(This,varargin)
            This = setoptions@report.tabularobj(This,varargin{:});
            This.options.long = false;
        end
        
    end

    methods (Access=protected,Hidden)
        
        varargout = mycompilepdf(varargin)
        varargout = mysubplot(varargin)
        varargout = myplot(varargin)
        varargout = speclatexcode(varargin)
        
    end
    
end