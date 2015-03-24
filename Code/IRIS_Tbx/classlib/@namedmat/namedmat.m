classdef namedmat < double
    % namedmat  Matrices with named rows and columns.
    %
    % Matrices with named rows and columns are returned by several IRIS
    % functions, such as [model/acf](model/acf), [model/xsf](model/xsf),
    % or [model/fmse](model/fmse), to facilitate easy selection of
    % submatrices by referrring to variable names in rows and columns.
    %
    % Namedmat methods:
    %
    % Constructor
    % ============
    %
    % * [`namedmat`](namedmat/namedmat) - Create a new matrix with named rows and columns.
    %
    % Manipulating named matrices
    % ============================
    %
    % * [`select`](namedmat/select) - Select submatrix by referring to row names and column names.
    % * [`transpose`](namedmat/transpose) - Transpose each page of matrix with names rows and columns.
    %
    % Getting row and column names
    % =============================
    %
    % * [`rownames`](namedmat/rownames) - Names of rows in namedmat object.
    % * [`colnames`](namedmat/colnames) - Names of columns in namedmat object.
    %
    % Sample characteristics
    % =======================
    %
    % * [`cutoff`](namedmat/cutoff] - 
    %
    % All operators and functions available for standard Matlab matrices
    % and arrays (i.e. double objects) are also available for namedmat
    % objects.
    %
    
    % -IRIS Toolbox.
    % -Copyright (c) 2007-2013 IRIS Solutions Team.
    
    properties (SetAccess = protected)
        Rownames = {};
        Colnames = {};
    end
    
    methods
        function this = namedmat(x,varargin)
            % namedmat  Create a new matrix with named rows and columns.
            %
            % Syntax
            % =======
            %
            %     X = namedmat(X,ROWNAMES,COLNAMES)
            %
            % Input arguments
            % ================
            %
            % * `X` [ numeric ] - Matrix or multidimensional array.
            %
            % * `ROWNAMES` [ cellstr ] - Names for individual rows of `X`.
            %
            % * `COLNAMES` [ cellstr ] - Names for individual columns of
            % `X`.
            %
            % Output arguments
            % =================
            %
            % * `X` [ namedmat ] - Matrix with named rows and columns.
            %
            % Description
            % ============
            %
            % Namedmat objects are used by some of the IRIS functions to
            % preserve the names of variables that relate to individual
            % rows and columns, such as in
            %
            % * `acf`, the autocovariance and autocorrelation functions,
            % * `xsf`, the power spectrum and spectral density functions,
            % * `fmse`, the forecast mean square error fuctions,
            % * etc.
            %
            % You can use the function [`select`](namedmat/select) to
            % extract submatrices by referring to a selection of names.
            %
            % Namedmat matrices derives from the built-in double class of
            % objects, and hence you can use any operators and functions on
            % them that are available for double objects.
            %
            % Example
            % ========
            %
            
            % -IRIS Toolbox.
            % -Copyright (c) 2007-2013 IRIS Solutions Team.
            
            %**************************************************************************
            
            if nargin == 0
                x = [];
            end
            this = this@double(x);
            if ~isempty(varargin)
                this.Rownames = varargin{1};
                varargin(1) = [];
                if length(this.Rownames) ~= size(x,1)
                    utils.error('namedmat', ...
                        'Number of rows must match number of row names.');
                end
            end
            if ~isempty(varargin)
                this.Colnames = varargin{1};
                varargin(1) = []; %#ok<NASGU>
                if length(this.Rownames) ~= size(x,1)
                    utils.error('namedmat', ...
                        'Number of columns must match number of column names.');
                end
            end
        end
        
        function disp(this)
            disp(double(this));
            addspace = false;
            if ~isempty(this.Rownames)
                disp(['   Rows:',sprintf(' %s',this.Rownames{:})]);
                addspace = true;
            end
            if ~isempty(this.Colnames)
                disp(['Columns:',sprintf(' %s',this.Colnames{:})]);
                addspace = true;
            end
            if addspace
                strfun.loosespace();
            end
        end

        varargout = colnames(varargin)
        varargout = cutoff(varargin);
        varargout = horzcat(varargin)
        varargout = plot(varargin)
        varargout = rownames(varargin)
        varargout = select(varargin)
        varargout = subsasgn(varargin)
        varargout = subsref(varargin)
        varargout = transpose(varargin)
        varargout = vertcat(varargin)
        
    end
    
end