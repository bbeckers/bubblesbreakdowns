% publish  Compile PDF from report object.
%
% Syntax
% =======
%
%     [OutpFile,Rerun] = P.publish(InpFile,...)
%
% Input arguments
% ================
%
% * `P` [ struct ] - Report object created by the `report.new` function.
%
% * `InpFile` [ char ] - File name under which the compiled PDF will be
% saved.
%
% Output arguments
% =================
%
% * `OutpFile` [ char ] - Name of the resulting PDF.
%
% * `Rerun` [ numeric ] -  Number of times the LaTeX compiler was run to
% compile the PDF.
%
% Options
% ========
%
% * `'abstract='` [ char | *empty* ] - Abstract that will displayed on the
% title page.
%
% * `'abstractWidth='` [ numeric | *`1`* ] - Width of the abstract on the
% page as a percentage of the full default width (between `0` and `1`).
%
% * `'author='` [ char | *empty* ] - List of authors on the title page
% separated with `\and` or `\\`.
%
% * `'cd='` [ `true` | *`false`* ] - If `true` do not use the `pdflatex`
% option `-include-directory` and instead change the directory (`cd`)
% temporarily to the location of the input file; this is a workaround for
% systems where `pdflatex` does not support the option
% `-include-directory`.
%
% * `'cleanup='` [ *`true`* | `false` ] - Delete all temporary files
% created when compiling the report.
%
% * `'compile='` [ *`true`* | `false` ] - Compile the source files to an
% actual PDF; if `false` only the source files are created.
%
% * `'date='` [ char | *`'\today'`* ] - Date on the title page.
%
% * `'display='` [ *`true`* | `false` ] - Display the \LaTeX compiler
% report on the final iteration.
%
% * `'echo='` [ `true` | *`false`* ] - If `true`, the optional flag
% `'-echo'` will be used in the Matlab function `system` when compiling the
% PDF; this causes the screen output and all prompts to be displayed
% for each run of the compiler.
%
% * `'epsToPdf='` [ char | *`Inf`* ] - Command line arguments for EPSTOPDF;
% `Inf` means OS-specific arguments are used.
%
% * `'fontEnc='` [ char | *`'T1'`* ] - \LaTeX\ font encoding.
%
% * `'makeTitle='` [ `true` | *`false`* ] - Produce title page (with title,
% author, date, and abstract).
%
% * `'package='` [ char | cellstr | *empty* ] - Package or list of packages
% that will be imported in the preamble of the LaTeX file.
%
% * `'paperSize='` [ `'a4paper'` | *`'letterpaper'`* ] - Paper size.
%
% * `'orientation='` [ *`'landscape'`* | `'portrait'` ] - Paper orientation.
%
% * `'preamble='` [ char | *empty* ] - \LaTeX\ commands that will be placed
% in the \LaTeX\ file preamble.
% 
% * `'timeStamp='` [ char | *`'datestr(now())'`* ] - String printed in the
% top-left corner of each page.
%
% * `'tempDir='` [ char | *`tempname(cd())`* ] - Directory for
% storing temporary files; the directory is deleted at the end of the
% execution if it's empty.
%
% * `'maxRerun='` [ numeric | *`5`* ] - Maximum number of times the \LaTeX\
% compiler will be run to resolve cross-references, etc.
%
% * `'minRerun='` [ numeric | *`1`* ] - Minimum number of times the \LaTeX\
% compiler will be run to resolve cross-references, etc.
%
% * `'textScale='` [ numeric | *`0.8`* ] - Percentage of the total page
% area that will be used.
%
% Description
% ============
% 
% Difference between `'display='` and `'echo='`
% ----------------------------------------------
%
% There are two differences between these otherwise similar options:
%
% * When publishing the final PDF, the PDFLaTeX compiler may be called more
% than once to resolve cross-references, the table of contents, and so on.
% Setting `'display=' true` only displays the screen output from the final
% iteration only, while `'echo=' true` displays the screen outputs from all
% iterations.
%
% * In the case of a compiler error unrelated to the \LaTeX\ code, the
% compiler may stop and prompt the user to respond. The prompt only appears
% on the screen when `'echo=' true`. Otherwise, Matlab may remain in a busy
% state with no on-screen information, and `Ctrl+C` may be needed to regain
% control.
%
% Example
% ========
%

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.
