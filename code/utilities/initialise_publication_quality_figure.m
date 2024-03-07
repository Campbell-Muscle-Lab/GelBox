function subplots=initialise_publication_quality_figure(varargin)
% Function maximizes figure size, sets up the print options

% Default values
params.figure_handle=gcf;
params.left_margin=0.5;             % 0.5 inch margin
params.right_margin=4.6;            % One column (3.5 inch) wide
params.top_margin=0.0;              % 0.5 inch margin
params.bottom_margin=0.5;
params.no_of_panels_wide=2;
params.no_of_panels_high=2;
params.axes_padding_left=0.75;  	% 0.75 inch left padding for labels
params.axes_padding_right=0.25;     % 0.25 inch right padding
params.axes_padding_top=0.25;		% 0.25 inch top padding
params.x_to_y_axes_ratio=1.0;       % ratio of x to y axes lengths
params.panel_label_font_size=14;    % font size label
params.font_name='Helvetica';
params.starting_letter=0;
params.individual_panel_labels='';
params.axes_padding_bottom=0.75;
params.relative_row_heights=[];
params.right_dead_space=0;
params.individual_padding=0;
params.left_pads=[];
params.right_pads=[];

params.left_subplot_adjustments=[];
params.right_subplot_adjustments=[];
params.bottom_subplot_adjustments=[];
params.height_subplot_adjustments=[];

params.panel_label_x_offset=[];
params.panel_label_y_offset=[];

params.individual_panels_wide=0;
params.omit_panels=[];


% Check for overrides
params=parse_pv_pairs(params,varargin);

% Updates
if (length(params.axes_padding_bottom)==1)
    params.axes_padding_bottom=params.axes_padding_bottom * ...
        ones(params.no_of_panels_high,1);
end

if (isempty(params.relative_row_heights))
    params.relative_row_heights=ones(params.no_of_panels_high,1);
end

if (params.individual_panels_wide == 0)
    params.individual_panels_wide = params.no_of_panels_wide * ...
                                        ones(1, params.no_of_panels_high);
end

% Error checking
if (length(params.axes_padding_bottom)~=params.no_of_panels_high)
    disp('Axes padding problem');
end

% Do some preparatory calculations
for row_counter=1:params.no_of_panels_high
    
    across(row_counter)=params.no_of_panels_wide;
    if (params.individual_panels_wide(1)>0)
        across(row_counter)=params.individual_panels_wide(row_counter);
    end
    
    if (~params.individual_padding)
        axes_width_inches(row_counter)= ...
            (8.5-(params.left_margin+params.right_margin+ ...
                params.right_dead_space)- ...
            (params.no_of_panels_wide*(params.axes_padding_left+ ...
                					params.axes_padding_right)))/ ...
            across(row_counter);
    else
        % Calculate the row indices
        if (row_counter==1)
            row_indices=1:across;
        else
            row_indices= ...
                sum(params.individual_panels_wide(1:(row_counter-1))) + ...
                (1:across(row_counter));
        end
        
        %  Override if horizontal paddings is individually specified
        axes_width_inches(row_counter) = (8.5 -  ...
            (params.left_margin+params.right_margin) - ...
            sum(params.left_pads(row_indices)) - ...
            sum(params.right_pads(row_indices))) / ...
                across(row_counter);
    end
end
        
                            
base_axes_height_inches=axes_width_inches(1)/params.x_to_y_axes_ratio;

% Figure width and height in inches
figure_width=8.5-(params.left_margin+params.right_margin);

if (length(params.axes_padding_top)==1)
    
    figure_height= params.top_margin + ...
        (params.no_of_panels_high * params.axes_padding_top) + ...
            sum(params.axes_padding_bottom(1:params.no_of_panels_high)) + ...
            ((sum(params.relative_row_heights(1:params.no_of_panels_high))) ...
                *base_axes_height_inches+ ...
        params.bottom_margin);
else
    figure_height = params.top_margin + ...
        (sum(params.axes_padding_top(1:params.no_of_panels_high))) + ...
            sum(params.axes_padding_bottom(1:params.no_of_panels_high)) + ...
            ((sum(params.relative_row_heights(1:params.no_of_panels_high))) ...
                *base_axes_height_inches+ ...
        params.bottom_margin);
end

ptm=params.top_margin;
pnph=params.no_of_panels_high;
paxpt=params.axes_padding_top;
papb=params.axes_padding_bottom;
prrh=params.relative_row_heights;
pbm=params.bottom_margin;

if (length(params.left_subplot_adjustments)==0)
    lhs_adjustments=zeros( ...
        params.no_of_panels_wide*params.no_of_panels_high,1);
else
    lhs_adjustments=params.left_subplot_adjustments;
end

if (length(params.right_subplot_adjustments)==0)
    rhs_adjustments=zeros( ...
        params.no_of_panels_wide*params.no_of_panels_high,1);
else
    rhs_adjustments=params.right_subplot_adjustments;
end


if (length(params.bottom_subplot_adjustments)==0)
    bot_adjustments=zeros( ...
        params.no_of_panels_wide*params.no_of_panels_high,1);
else
    bot_adjustments=params.bottom_subplot_adjustments;
end

if (length(params.height_subplot_adjustments)==0)
    height_adjustments=zeros( ...
        params.no_of_panels_wide*params.no_of_panels_high,1);
else
    height_adjustments=params.height_subplot_adjustments;
end

% Clear figure
figure(params.figure_handle);
clf;
set(params.figure_handle,'Units','inches','PaperType','usletter');
set(params.figure_handle,'Position', ...
    [params.left_margin 9-figure_height figure_width figure_height]);

% Loop through sub-plots

subplots=[];
subplot_counter=0;

for row_counter=1:params.no_of_panels_high
  
    % Calculate the row indices
    if (row_counter==1)
        row_indices=1:across;
    else
        if (params.individual_panels_wide(1)>0)
            row_indices= ...
                sum(params.individual_panels_wide(1:(row_counter-1))) + ...
                (1:across(row_counter));
        else
            row_indices=(row_counter-1)*params.no_of_panels_wide + ...
                (1:params.no_of_panels_wide);
        end
    end
   
    if (~params.individual_panels_wide(1))
        holder=params.no_of_panels_wide;
    else
         holder=params.individual_panels_wide(row_counter);
    end
        
    for column_counter=1:holder
        
        subplot_counter=subplot_counter+1;
        
        % Check for omit panels
        if (any(params.omit_panels==subplot_counter))
            continue;
        end
        
        % Set lhs normalized to figure width
        if (~params.individual_padding)
            lhs=((column_counter-1)* ...
                (params.axes_padding_left+axes_width_inches(row_counter)+ ...
                        params.axes_padding_right) + ...
                params.axes_padding_left)/figure_width;
        else
            lhs = (((column_counter-1)*axes_width_inches(row_counter) ) + ...
                sum(params.left_pads(row_indices(1:column_counter))) + ...
                sum(params.right_pads(row_indices(1:column_counter-1)))) / figure_width;
        end
   
        if (length(params.axes_padding_top)==1)
            bottom=(figure_height - params.top_margin - ...
                        (row_counter*params.axes_padding_top) -...
                        (sum(params.relative_row_heights(1:row_counter))* ...
                            base_axes_height_inches) - ...
                        bot_adjustments(subplot_counter) - ...
                        sum(params.axes_padding_bottom(1:row_counter-1)))/ ...
                    figure_height;
        else
            bottom=(figure_height - params.top_margin - ...
                        (sum(params.axes_padding_top(1:row_counter))) -...
                        (sum(params.relative_row_heights(1:row_counter))* ...
                            base_axes_height_inches - ...
                            bot_adjustments(subplot_counter)) - ...
                        sum(params.axes_padding_bottom(1:row_counter-1)))/ ...
                    figure_height;
        end
            
        subplots(subplot_counter)=subplot('Position', ...
            [lhs+(lhs_adjustments(subplot_counter)/figure_width) bottom ...
            (axes_width_inches(row_counter)- ...
                lhs_adjustments(subplot_counter)- ...
                rhs_adjustments(subplot_counter))/figure_width ...
                (params.relative_row_heights(row_counter)* ...
                    base_axes_height_inches + ...
                    bot_adjustments(subplot_counter) + ...
                    height_adjustments(subplot_counter))/figure_height]);
                
        % Draw Label
        if (params.panel_label_font_size>0)
            
            if (0)
            % Create the subplot and set up the coordinates for the label
            subplot(subplots(subplot_counter));
                        
            x_pos=-params.axes_padding_left;
            y_pos=params.axes_padding_top + ...
                params.relative_row_heights(row_counter)* ...
                    base_axes_height_inches;
            if (row_counter==1)
                y_pos=y_pos+params.top_margin;
            end
            text(x_pos,y_pos, ...
                char(subplot_counter+64+params.starting_letter), ...
                'FontSize',params.panel_label_font_size, ...
                'FontWeight','bold', ...
                'Units','inches', ...
                'HorizontalAlignment','left', ...
                'VerticalAlignment','top', ...
                'FontName',params.font_name ...
                );
            
            text('Units','data');
            else
                % Move to the subplot
                h=subplot(subplots(subplot_counter));
                set(h,'Units','inches');
                pos_vector=get(h,'Position');
                lhs=pos_vector(1);
                top=pos_vector(2)+pos_vector(4);
                
                % Find the positions for the labels
                if (~params.individual_padding)
                    x_pos=-params.axes_padding_left;
                else
                    x_pos=-params.left_pads( ...
                        row_indices(column_counter));
                end
                x_pos=x_pos-lhs_adjustments(subplot_counter);
                
                if (~isempty(params.panel_label_x_offset))
                    x_pos = x_pos + params.panel_label_x_offset(subplot_counter);
                end
                
                if (length(params.axes_padding_top)==1)
                    y_pos=params.axes_padding_top + ...
                        params.relative_row_heights(row_counter)* ...
                            base_axes_height_inches;
                else
                    y_pos=params.axes_padding_top(row_counter) + ...
                        params.relative_row_heights(row_counter)* ...
                            base_axes_height_inches;
                end
                if (row_counter==1)
                    y_pos=y_pos+params.top_margin;
                end
                y_pos=y_pos+bot_adjustments(subplot_counter)+ ...
                    height_adjustments(subplot_counter);
                
                if (~isempty(params.panel_label_y_offset))
                    y_pos = y_pos + params.panel_label_y_offset(subplot_counter);
                end
                
                                
                if (length(params.individual_panel_labels)>0)
                    text_string=params.individual_panel_labels{ ...
                        subplot_counter};
                else
                    text_string=sprintf('%c', ...
                        subplot_counter+64+params.starting_letter);
                end
                
                % Slight x offset added to prevent labels being cropped
                text(x_pos+0.01,y_pos,text_string, ...
                    'Units','inches', ...
                    'FontSize',params.panel_label_font_size, ...
                    'FontWeight','bold', ...
                    'Units','inches', ...
                    'HorizontalAlignment','left', ...
                    'VerticalAlignment','top', ...
                    'FontName',params.font_name ...
                    );
                text('Units','data');
            end
        end
    
        hold on;
        drawnow;
    end
end

% Restore defaults
set(gcf,'PaperUnits','inches')

        
