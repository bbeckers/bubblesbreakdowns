% !transition_equations  Block of transition equations.
%
% Syntax
% =======
%
%     !transition_equations
%         Equation1;
%         Equation2;
%         Equation2;
%         ...
% Short-cut syntax
% =================
%
%     !equations
%         Equation1;
%         Equation2;
%         Equation3;
%         ...
%
% Syntax with equation labels
% ============================
%
%     !transition_equations
%         Equation1;
%         'Equation label' Equation2;
%         Equation3;
%         ...
%
% Description
% ============
%
% The `!transition_equations` keyword starts a new block of transition
% equations (i.e. endogenous equations); the eqautions can stretch over
% multiple lines and must be separated by semi-colons. You can have as many
% equation blocks as you wish in any order in your model file: They all
% get combined together when you read the model file in.
% 
% You can add descriptive labels to the equations (in single or double
% quotes, preceding the equation); these will be stored in, and
% accessible from, the model object.
%
% Example
% ========
%
%     !transition_equations
%         'Euler equation' C{1}/C = R*beta;

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.
