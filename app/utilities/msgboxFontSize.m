function [th, bth] = msgboxFontSize(varargin)
% [th, bth] = msgboxFontSize(h, fontsize, varargin)
% The text properties in Matlab's msgbox(), errordlg(), warndlg(), etc cannot be directly 
% edited and the default fontsize is quite small.  This function gets around that problem 
% by searching for text within the msgbox handle, changes its fontsize, and then changes 
% the window size to fit the new fontsize.  The position of the lower, left corner of 
% the msgbox will not change and the resize property is turned on.  Additional name-value
% text properties can be set except for 'units' (see examples).  
%
% To access the icons in errordlg, warndlg, etc, if h is the handle to the dialog box, 
%   iconAxes = findall(h,'tag','IconAxes'); 
%   iconAxes.Position([3,4]) = iconAxes.Position([3,4]) + [6,6]; %increase size
%
% This function has been tested in r2014a, r2016a, r2017b, r2019a. Please direct all 
% problems and high fives to the email address at the bottom of this text.
%
% INPUT
%   h = handle(s) to dialog that is already created. 
%   fonsize = (optional) the fontsize you'd like. 
%             Default is 12; FONTSIZE should either be a 
%             scalar value to be applied to all handles
%             in h or an array the same size of h to 
%             specify the fontsize for each element in h.
% NAME-VALUE PAIRS
%  * 'ignorebuttons': true/false or 1/0; when true, the changes
%       will not be made to the buttons within the dialog. Note
%       that button size does not change (for now) so large 
%       fonts may exceed the button space.  
%  * Any text name-value pair can be included (see link below).
%       https://www.mathworks.com/help/matlab/ref/matlab.graphics.primitive.text-properties.html
% OUTPUT
%   th = handles to the text objects in same order as 'h'.
%   bth = a vector of button handles so the user can manually
%       edit their properties.
%
% EXAMPLES
%   h = msgbox('Hello world (20 pt font)', 'Example');
%   msgboxFontSize(h,20); 
%   msgboxFontSize(h,20,'ignorebuttons',true);
%
%   str = 'For longer messages, the window will expand to fit the text!'; 
%   h = msgbox(str, 'Longer example');
%   msgboxFontSize(h, 18); 
%
%   str = sprintf('List 1\nList 2\nList 3\nList 4\nList 5\nList 5'); 
%   h = msgbox(str, 'Tall example'); 
%   msgboxFontSize(h, 18)
%
%   h1 = msgbox('box 1', 'Small font'); 
%   h2 = msgbox('box 2', 'Large font'); 
%   msgboxFontSize([h1, h2], 14); 
%   msgboxFontSize([h1, h2], [10, 16]); 
%
%   edh = errordlg('Error dialog message', 'Error name'); 
%   msgboxFontSize(edh,12,'FontName','Consolas','Color','r','FontSize',14); 
% 
% Danz 180810
% Copyright (c) 2018, Adam Danz 
% All rights reserved.
% adam.danz@gmail.com

% Updates
% 190625 added ability to input name-value text params;
%        and added input parser. (v 2.0)
% 190708 prevent set() diplay; added button handle output.  (v 2.1)
% 190724 updated to work in r2014a (and possibly earlier releases)

%% Parse inputs
p = inputParser(); 
p.FunctionName = mfilename;
p.KeepUnmatched = true;         %accept additional parameter value inputs 

checkIgnoreButtons = @(x)(islogical(x) || isnumeric(x)); 
checkhandles = @(h)all(ishghandle(h));
addRequired(p, 'h', checkhandles);
addOptional(p, 'fontsize', 12, @isnumeric);
addParameter(p, 'ignorebuttons', false, checkIgnoreButtons);
parse(p,varargin{:})

% Prepare the unmatched text() parameters.
% If a param is passed that isn't accepted by text(), an error is thrown from text() function.
unmatchNameVal = reshape([fieldnames(p.Unmatched)'; struct2cell(p.Unmatched)'], 1, []); 

%% Input validation
% replicate fontsize if needed
if numel(p.Results.h) ~= numel(p.Results.fontsize)
    if numel(p.Results.fontsize) == 1
        fontsize = repmat(p.Results.fontsize, size(p.Results.h)); 
    else
        error('FONTSIZE input must either be a scalar or an array the size of h')
    end
else
    fontsize = p.Results.fontsize; 
end

% If user specified units, remove and throw warning.
unitsCheck = strcmpi(unmatchNameVal,'Units'); 
if any(unitsCheck)
    unmatchNameVal(find(unitsCheck)+[0,1]) = []; 
    warning('The "Units" name-value pair was ignored.')
end

% If user wants to apply changes to buttons, convert 'color' to 'ForegroundColor'
if ~p.Results.ignorebuttons
    unmatchNameValBtn = regexprep(unmatchNameVal,'Color','ForegroundColor','ignorecase'); 
end
  
%% Adjust fontsize and window size
minSize = [150, 51]; %units = points; %the smallest allowed window
th = gobjects(size(p.Results.h)); 
for i = 1:length(p.Results.h)
    origUnits = get(p.Results.h(i),'Units'); 
    set(p.Results.h(i),'Units','points'); 
    % Find text object
    th(i) = findall(p.Results.h(i), 'Type', 'Text');
    % set optional name-value text properties
    if ~isempty(unmatchNameVal)
        set(th(i), unmatchNameVal{:});
    end 
    set(th(i),'FontSize',fontsize(i)); % see footnote [2]
    thExt = get(th(i),'Extent');
    hPos = get(p.Results.h(i),'Position'); 
    deltaWidth = sum(thExt([1,3]))-hPos(3) + thExt(1);
    deltaHeight = sum(thExt([2,4]))-hPos(4) + 10;
    set(p.Results.h(i),'Position',[hPos(1:2), max(hPos([3,4]) + [deltaWidth, deltaHeight], minSize)]); 
    % apply to buttons unless ignored
    if ~p.Results.ignorebuttons
        btnHand = findall(p.Results.h(i), 'Style', 'pushbutton');
        if ~isempty(unmatchNameVal)
            set(btnHand, unmatchNameValBtn{:});
        end
        set(btnHand,'FontSize',fontsize(i));
    end
    % Allow user to resize
    set(p.Results.h(i),'Resize','on'); % See footnote [1]
    set(p.Results.h(i),'Units',origUnits); 
end

%% Output handles to buttons (if any) so user can manually alter properties.
bth = findall(p.Results.h, 'Style', 'pushbutton'); 

%% Footnotes
% [1] in r2014a and possibly other older releases, setting 'Resize' 'on'
%   additionally increases the figure size a bit.  Why? Who knows. 
% [2] in r2014a and possibly other older releases, increasing the fontsize
%   results in a pixelated but readable font. 
