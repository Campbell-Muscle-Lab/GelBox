function axes_data=improve_axes(varargin)
% Makes high quality axes by deleting defaults and redrawing from scratch

% Set defaults
params.axis_handle=[];
params.x_ticks=[];
params.y_ticks=[];
params.axis_line_width=1.5;
params.x_axis_label='X axis';
params.y_axis_label='Y axis';
params.label_font_size=12;
params.font_name='Arial';
params.x_limits=[];
params.x_tick_labels=[];
params.omit_x_tick_labels=0;
params.x_tick_label_positions=[];
params.x_tick_label_rotation=0;
params.y_tick_label_rotation=0;
params.y_tick_label_positions=[];
params.y_tick_labels={};
params.x_tick_decimal_places=3;
params.y_tick_decimal_places=3;
params.tick_font_size=12;
params.x_tick_label_font_size = [];
params.x_log_mode=0;
params.y_log_mode=0;
params.x_axis_colour=[0 0 0];
params.x_axis_off=0;
params.y_axis_off=0;
params.x_tick_rotation=0;
params.redraw_axes=1;
params.x_axis_offset=0.05;
params.x_label_offset=-0.25;
params.y_axis_offset=0.025;
params.y_label_offset=-0.25;
params.y_label_vertical_offset=0;
params.axis_color='k';
params.y_tick_label_horizontal_offset=-0.04;
params.y_tick_length=0.025;
params.x_tick_label_vertical_offset=-0.05;
params.x_tick_label_horizontal_offset=0;
params.x_tick_length=[];
params.x_ticks_off=0;
params.x_tick_dir=1;
params.y_tick_dir=1;
params.y_label_rotation=0;
params.x_label_rotation=0;
params.y_label_text_interpreter='tex';
params.x_label_text_interpreter='tex';
params.x_tick_label_text_interpreter='tex';
params.title='';
params.title_x_offset=NaN; % NaN is center, number defines placement
params.title_y_offset=1.05;
params.title_v_Align='middle';
params.title_h_Align='center';
params.title_font_size=12;
params.title_font_weight='normal';
params.title_font_angle='normal';
params.title_text_interpreter='tex';
params.title_edge_color='none';
params.panel_label='';
params.panel_label_x_offset=0;
params.panel_label_y_offset=0;
params.panel_label_font_size=12;
params.y_adjust = 1;                % scales to match x axes
params.gui_scale_factor = 0;
params.clip_x_axis = 0;

% Update
params=parse_pv_pairs(params,varargin);

if (isempty(params.axis_handle))
    params.axis_handle = gca;
end

% Wipe out the old axis
set(params.axis_handle,'Visible','off');

% Set the limits if they are not specified
if (isempty(params.x_ticks))
    params.x_ticks=xlim(params.axis_handle);
end
if (isempty(params.y_ticks))
    params.y_ticks=ylim(params.axis_handle);
end

% Try to set tick decimal places intelligently
if (isempty(params.y_tick_decimal_places))
    n=log10(params.y_ticks(end));
    if (n<=0)
        params.y_tick_decimal_places = -round(n)+2;
    else
        if (n>=2)
            params.y_tick_decimal_places = 0;
        else
            params.y_tick_decimal_places = 1;
        end
    end
end

% Let's do easy things first
y_tick_length = params.y_tick_length * ...
    abs((params.x_ticks(end)-params.x_ticks(1)));
% Set y axis position
y_axis_x_location = params.x_ticks(1)- ...
    params.y_axis_offset*abs((params.x_ticks(end)-params.x_ticks(1)));

% Set figure x range
if (params.x_ticks(end)>params.x_ticks(1))
    lhs=y_axis_x_location - y_tick_length;
else
    lhs=y_axis_x_location + y_tick_length;
end
rhs=params.x_ticks(end);

if (params.y_ticks(end)>params.y_ticks(1))
    % Increasing y axis

    % First work out the x tick length
    if (~isempty(params.x_tick_length))
        x_tick_length = params.x_tick_length * ...
            (params.y_ticks(end)-params.y_ticks(1));
    else
        % Set x_tick length to the same relative length as the y ticks
        rel_y=params.y_tick_length;
        pos=get(params.axis_handle,'Position');
        
        x_tick_length = (params.y_ticks(end)-params.y_ticks(1))* ...
            params.y_adjust * ...
            rel_y*(pos(3)/pos(4));
    end
   
    % Set x axis position
    x_axis_y_location = params.y_ticks(1) - ...
        params.x_axis_offset*(params.y_ticks(end)-params.y_ticks(1));
    
    % Set figure y range
    bottom = x_axis_y_location - x_tick_length;
    top = params.y_ticks(end);
else
    % Decreasing y axis
    
    % First work out the x tick length
    if (~isempty(params.x_tick_length))
        x_tick_length = params.x_tick_length * ...
            (params.y_ticks(end)-params.y_ticks(1));
    else
        % Set x_tick length to the same relative length as the y ticks
        rel_y=params.y_tick_length;
        pos=get(params.axis_handle,'Position');
        
        x_tick_length = (params.y_ticks(end)-params.y_ticks(1))* ...
            rel_y*(pos(3)/pos(4));
    end
    
    % Set x axis position
    x_axis_y_location = params.y_ticks(1) - ...
        params.x_axis_offset*(params.y_ticks(end)-params.y_ticks(1));
    
    % Set figure y range
    bottom = x_axis_y_location - x_tick_length;
    top = params.y_ticks(end);
end

% Draw y axis
if (~params.y_axis_off)
    line([lhs y_axis_x_location y_axis_x_location lhs], ...
        [params.y_ticks(1)*[1 1] params.y_ticks(end)*[1 1]], ...
        'LineWidth',params.axis_line_width, ...
        'Color',params.axis_color, ...
        'Clipping','off', ...
        'Parent',params.axis_handle);

    if (isempty(params.y_tick_labels))
        % Draw the additional default y_axis ticks
        if (length(params.y_ticks)>2)
            for tick_counter=2:length(params.y_ticks)-1
                line([lhs y_axis_x_location], ...
                    params.y_ticks(tick_counter)*[1 1], ...
                    'LineWidth',params.axis_line_width, ...
                    'Color',params.axis_color, ...
                    'Parent',params.axis_handle);
            end
        end
    else
        % Draw ticks at y_tick positions
        for tick_counter=1:length(params.y_tick_label_positions)
            y_tick_label_position= ...
                params.y_tick_label_positions(tick_counter);
            if ((y_tick_label_position~=params.y_ticks(1))&& ...
                    (y_tick_label_position~=params.y_ticks(end)))
                line([lhs y_axis_x_location], ...
                    y_tick_label_position*[1 1], ...
                    'LineWidth',params.axis_line_width, ...
                    'Color',params.axis_color, ...
                    'Parent',params.axis_handle);
            end
        end
    end
end

% Draw x axis
if (~params.x_axis_off)
    if(params.clip_x_axis)
        c = 'w';
    else
        c = params.axis_color;
    end
    line([params.x_ticks(1)*[1 1] params.x_ticks(end)*[1 1]], ...
        [bottom x_axis_y_location x_axis_y_location bottom], ...
        'LineWidth',params.axis_line_width, ...
        'Color',c, ...
        'Clipping','off', ...
        'Parent',params.axis_handle);
    
    if(params.clip_x_axis)
        
        line([params.x_tick_label_positions(1)*[1 1] params.x_tick_label_positions(end)*[1 1]], ...
            [bottom x_axis_y_location x_axis_y_location bottom], ...
            'LineWidth',params.axis_line_width, ...
            'Color',params.axis_color, ...
            'Clipping','off', ...
            'Parent',params.axis_handle);
    end
        
        
    

    % Draw the additional x_axis ticks
    if (isempty(params.x_tick_labels))
        % Default behavior
        if (length(params.x_ticks)>2)
            for tick_counter=2:length(params.x_ticks)-1
                line(params.x_ticks(tick_counter)*[1 1], ...
                    [bottom x_axis_y_location], ...
                    'LineWidth',params.axis_line_width, ...
                    'Color',params.axis_color, ...
                    'Parent',params.axis_handle);
            end
        end
    else
        % Draw ticks at z_tick positions
        for tick_counter=1:length(params.x_tick_label_positions)
            x_tick_label_position= ...
                params.x_tick_label_positions(tick_counter);
            if ((x_tick_label_position~=params.x_ticks(1))&& ...
                    (x_tick_label_position~=params.x_ticks(end)))
                line(x_tick_label_position*[1 1], ...
                    [bottom x_axis_y_location], ...
                    'LineWidth',params.axis_line_width, ...
                    'Color',params.axis_color, ...
                    'Parent',params.axis_handle);
            end
        end
    end
end

% Draw y axis label
if (~params.y_axis_off)
    offset=params.y_label_offset*(params.x_ticks(end)-params.x_ticks(1));
    
    y_axis_label_x_location = y_axis_x_location + offset;
    
    text(y_axis_label_x_location, ...
        mean([params.y_ticks(1) params.y_ticks(end)]) + ...
            params.y_label_vertical_offset* ...
                (params.y_ticks(end)-params.y_ticks(1)), ...
        params.y_axis_label, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontName',params.font_name, ...
        'FontSize',params.label_font_size, ...
        'Color',params.axis_color, ...
        'Rotation',params.y_label_rotation, ...
        'Interpreter',params.y_label_text_interpreter, ...
        'Parent',params.axis_handle);
    
    

    % Draw y tick labels
    % Deduce the offset
    offset=params.y_tick_label_horizontal_offset* ...
        (params.x_ticks(end)-params.x_ticks(1));
     % Set alignment
    if (params.y_tick_label_rotation==0)
        hor_align='right';
        ver_align='middle';
    else
        hor_align='right';
        ver_align='top';
    end
    if (isempty(params.y_tick_labels))
        % Default tick labels
        for tick_counter=1:length(params.y_ticks)
            % Tick labels
            tick_string=print_number_to_specified_decimal_places( ...
               params.y_ticks(tick_counter),params.y_tick_decimal_places);

            % Draw
            text(y_axis_x_location+offset, ...
                params.y_ticks(tick_counter), ...
                tick_string, ...
                'FontName',params.font_name, ...
                'FontSize',params.tick_font_size, ...
                'Color',params.axis_color, ...
                'Rotation',params.y_tick_label_rotation, ...
                'HorizontalAlignment',hor_align, ...
                'VerticalAlignment',ver_align, ...
                'Parent',params.axis_handle);
        end
    else
        % User-provided ticks
        for tick_counter=1:length(params.y_tick_label_positions)
            % Draw
            text(y_axis_x_location+offset, ...
                params.y_tick_label_positions(tick_counter), ...
                params.y_tick_labels{tick_counter}, ...
                'FontName',params.font_name, ...
                'FontSize',params.tick_font_size, ...
                'Color',params.axis_color, ...
                'Rotation',params.y_tick_label_rotation, ...
                'HorizontalAlignment',hor_align, ...
                'VerticalAlignment',ver_align, ...
                'Parent',params.axis_handle);
        end
    end
end

% Draw x axis label
if (~params.x_axis_off)
    offset=params.x_label_offset*(params.y_ticks(end)-params.y_ticks(1));

    text((params.x_ticks(1)+params.x_ticks(end))/2, ...
        x_axis_y_location+offset, ...
        params.x_axis_label, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontName',params.font_name, ...
        'FontSize',params.label_font_size, ...
        'Color',params.axis_color, ...
        'Rotation',params.x_label_rotation, ...
        'Interpreter',params.x_label_text_interpreter, ...
        'Parent',params.axis_handle);

    % Draw x tick labels
    if (~params.omit_x_tick_labels)
        
        if (isempty(params.x_tick_label_font_size))
            font_size = params.tick_font_size;
        else
            font_size = params.x_tick_label_font_size;
        end
        
        % Deduce the offset
        offset=params.x_tick_label_vertical_offset* ...
            (params.y_ticks(end)-params.y_ticks(1));
         % Set alignment
        if (params.x_tick_label_rotation==0)
            hor_align='center';
            ver_align='top';
        else
            hor_align='right';
            ver_align='top';
        end
        if (isempty(params.x_tick_labels))
            % Default labels
            for tick_counter=1:length(params.x_ticks)
                tick_string=print_number_to_specified_decimal_places( ...
                   params.x_ticks(tick_counter),params.x_tick_decimal_places);

                % Draw
                text(params.x_ticks(tick_counter), ...
                    x_axis_y_location+offset, ...
                    tick_string, ...
                    'FontName',params.font_name, ...
                    'FontSize',font_size, ...
                    'Color',params.axis_color, ...
                    'Rotation',params.x_tick_label_rotation, ...
                    'HorizontalAlignment',hor_align, ...
                    'VerticalAlignment',ver_align, ...
                    'Parent',params.axis_handle);
            end
        else
            % User provided labels
            for tick_counter=1:length(params.x_tick_labels)
                % Draw
                text(params.x_tick_label_positions(tick_counter) + ...
                        params.x_tick_label_horizontal_offset, ...
                    x_axis_y_location+offset, ...
                    params.x_tick_labels(tick_counter), ...
                    'FontName',params.font_name, ...
                    'FontSize',font_size, ...
                    'Color',params.axis_color, ...
                    'Rotation',params.x_tick_label_rotation, ...
                    'HorizontalAlignment',hor_align, ...
                    'VerticalAlignment',ver_align, ...
                    'Interpreter',params.x_tick_label_text_interpreter, ...
                    'Parent',params.axis_handle);
            end
        end
    end
end

% Adjust scaling
f = params.gui_scale_factor;
new_lhs = lhs-f*(rhs-lhs);
new_rhs = rhs+f*(rhs-lhs);
new_bottom = bottom-f*(top-bottom);
new_top = top+f*(top-bottom);

if (new_top==new_bottom)
    new_top = new_bottom+1;
end

% Set the axes
if (new_rhs>new_lhs)
    xlim(params.axis_handle,[new_lhs new_rhs]);
else
    xlim(params.axis_handle,[new_rhs new_lhs]);
    set(params.axis_handle,'XDir','reverse');
end
if (new_bottom<new_top)
    ylim(params.axis_handle,[new_bottom new_top]);
else
    ylim(params.axis_handle,[new_top new_bottom]);
    set(params.axis_handle,'YDir','reverse');
end

% Draw title
% center title if title_x_offset=NaN

if (isnan(params.title_x_offset))
    text(mean([params.x_ticks(1) params.x_ticks(end)]), ...
        x_axis_y_location + (params.title_y_offset * diff(ylim(params.axis_handle))), ...
        params.title, ...
        'VerticalAlignment',params.title_v_Align, ...
        'HorizontalAlignment',params.title_h_Align, ...
        'FontName',params.font_name, ...
        'FontSize',params.title_font_size, ...
        'FontWeight',params.title_font_weight, ...
        'FontAngle',params.title_font_angle, ...
        'Color',params.axis_color, ...
        'Interpreter',params.title_text_interpreter, ...
        'EdgeColor',params.title_edge_color, ...
        'Parent',params.axis_handle);
else
    text( ...
        y_axis_x_location + (params.title_x_offset *diff(xlim(params.axis_handle))), ...
        x_axis_y_location + (params.title_y_offset *diff(ylim(params.axis_handle))), ...
        params.title, ...
        'VerticalAlignment',params.title_v_Align, ...
        'HorizontalAlignment',params.title_h_Align, ...
        'FontName',params.font_name, ...
        'FontSize',params.title_font_size, ...
        'FontWeight',params.title_font_weight, ...
        'FontAngle',params.title_font_angle, ...
        'Color',params.axis_color, ...
        'Interpreter',params.title_text_interpreter, ...
        'EdgeColor',params.title_edge_color, ...
        'Parent',params.axis_handle);
end


% Draw label

% Find the position of the subplot
h=subplot(params.axis_handle);
original_units=get(params.axis_handle,'Units');

set(h,'Units','inches');
pos_vector=get(h,'Position');
set(h,'Units',original_units);

lhs=pos_vector(1);
top=pos_vector(4);

x_pos=params.panel_label_x_offset;
y_pos=top+params.panel_label_y_offset;

% Now draw it
text(x_pos,y_pos,params.panel_label, ...
    'Units','inches', ...
    'FontSize',params.panel_label_font_size, ...
    'FontWeight','bold', ...
    'Units','inches', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','middle', ...
    'FontName',params.font_name, ...
    'Interpreter',params.title_text_interpreter, ...
    'Parent',params.axis_handle);
text('Units','data', ...
        'Parent',params.axis_handle);

% Set axes_data
axes_data.x_ticks = params.x_ticks;
axes_data.x_axis_y_location=x_axis_y_location;
axes_data.axis_line_width=params.axis_line_width;
axes_data.y_ticks = params.y_ticks;
axes_data.bottom = bottom;
if params.y_axis_off ~=1
axes_data.y_axis_label_x_location = y_axis_label_x_location;    
end
