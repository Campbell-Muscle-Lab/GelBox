function status=write_structure_to_excel(varargin)
% Written by Ken
% Code takes a structure and writes it to Excel
% NaNs in the original structure are written as strings
% Conventional strings are written correctly into the Excel file

params.filename = '';
params.structure = [];
params.sheet = 'Sheet1';

% Update
params=parse_pv_pairs(params,varargin);

% Code
status=0;

% Get full path
[path_string,name,ext]=fileparts(params.filename);
% If the path doesn't start with a drive, make it so
if ((numel(path_string)==0)||(~strcmp(path_string(2:3),':\')))
    params.filename = fullfile(cd,params.filename);
end

% Display
display_string = sprintf('Writing to %s: %s',params.filename,params.sheet);
disp(display_string);

% Suppress Excel warning
warning off MATLAB:xlswrite:AddSheet

% Create a cell array from the structure
cell_array = structure_to_cell_array(params.structure);

% Add field names to the structure
cell_array = [fieldnames(params.structure)' ; cell_array];

% First of all, check to see whether the file exists
if (exist(params.filename,'file'))
    [status,sheet_names]=xlsfinfo(params.filename);

    if (strcmp(status,'Microsoft Excel Spreadsheet'))
        % The file exists - we have to clear the appropriate sheet before
        % writing data to it.

        % We do this using activex

        Excel=actxserver('Excel.Application');
        % Make it invisible
        set(Excel,'Visible',0);
        % Turn off alerts
        set(Excel,'DisplayAlerts',0);
        % Get a handle to Excel's Workbooks
        Workbooks=Excel.Workbooks;
        % Open an Excel Workbook and activate it
        Workbook=Workbooks.Open(params.filename);
        all_sheets=Excel.ActiveWorkBook.Sheets;

        matching_index=find(strcmp(sheet_names,params.sheet));

        if (matching_index)
            current_sheet=get(all_sheets,'Item',matching_index);
            invoke(current_sheet.Rows,'Delete');
        end

        Workbook.Save;
        Workbooks.Close;
        invoke(Excel,'Quit');
    else
        % File isn't readable, delete it
        delete(params.filename);
    end
end
       
% Now write the new data
[status,message]=xlswrite(params.filename,cell_array,params.sheet);

% Error check
if (~status)
    msgbox(message.message,'Excel file problem','error');
end

end


function c = structure_to_cell_array(s)
% Function takes a structure and writes it to a cell array
% NaNs are written as 'NaN' strings

field_names=fieldnames(s);
no_of_field_names=length(field_names);
no_of_entries=0;
for field_counter=1:no_of_field_names
    no_of_field_entries(field_counter) = ...
        length(s.(field_names{field_counter}));
end

% Pre-allocate
c=cell(no_of_entries,no_of_field_names);

% Nested loop
for field_counter=1:no_of_field_names
    for entry_counter=1:no_of_field_entries(field_counter)
        x=s.(field_names{field_counter})(entry_counter);
        % Is it a cell
        if (iscell(x))
            if (isnumeric(cell2mat(x)))
                if (isnan(cell2mat(x)))
                    c{entry_counter,field_counter}='NaN';  
                else
                    c{entry_counter,field_counter}=x;
                end
            else
                % Check for leading digits in the string which create
                % issues in excel. If they are there, add in an extra '_'
                x=char(x);
%                 if (~isempty(x))
%                     if (length(str2num((x(1))))>0)
%                         x=['_' x];
%                     end
%                 end
                c{entry_counter,field_counter}=sprintf('%s',x);
            end
        else
            if (isnan(x))
                c{entry_counter,field_counter}='NaN';
            else
                c{entry_counter,field_counter}=x;
            end
        end
    end
            
end

end