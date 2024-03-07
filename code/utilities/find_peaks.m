function peak_data=find_peaks(varargin);

params.x=[];
params.y=[];
params.smoothing_half_window=0;
params.min_rel_delta_y=0.1;
params.min_x_index_spacing=20;

% Update values
params=parse_pv_pairs(params,varargin);

% Error checking
if (length(params.x)==0)
    error('No x data specified');
end
if (length(params.y)==0)
    error('No y data specified');
end
if (length(params.x)~=length(params.y))
    error('x and y data are unequal lengths');
end

% Calculations
min_y=min(params.y);
max_y=max(params.y);
min_real_delta_y=params.min_rel_delta_y*(max_y-min_y);
no_of_points=length(params.x);

% Smoothing if required
if (params.smoothing_half_window>0)
    params.y=smooth_array_in_place(params.y,params.smoothing_half_window);
end

% Start at the left-hand side
current_min_y=params.y(1);
current_min_index=1;
current_max_y=params.y(1);
current_max_index=1;
start_y=params.y(1);
start_position=1;

index=1;
last_extreme=0;
min_indices=[];
max_indices=[];


while (index<no_of_points)
    index=index+1;

    % Update the extreme holders
    if (params.y(index)>current_max_y)
        current_max_y=params.y(index);
        current_max_index=index;
    end
    if (params.y(index)<current_min_y)
        current_min_y=params.y(index);
        current_min_index=index;
    end

    % Have we gone far enough down and back up to mark a minimum?
    if ((params.y(index)-current_min_y)>min_real_delta_y)
        if ((current_max_y-current_min_y)>min_real_delta_y)
            % Yes - now, if this is the first extreme, check that you've
            % gone the requisite distance below the first point
            if ((last_extreme==0)& ...
                    ((params.y(1)-current_min_y)>min_real_delta_y)) | ...
                (last_extreme==1)
                % Are we far enough away from the minimum to know for sure
                if ((last_extreme==0) & ...
                            (index>params.min_x_index_spacing)) | ...
                    ((last_extreme==1) && ...
                            (index-current_min_index) > ...
                                params.min_x_index_spacing)
                
                    % We've defined a minimum
                    min_indices=[min_indices current_min_index];
                    current_min_y=params.y(index);
                    current_min_index=index;
                    current_max_y=params.y(index);
                    current_max_index=index;
                    last_extreme=-1;
                end
            end
        end
    end
    
    % Have we gone far enough up and back down to mark a maximum?
    if ((current_max_y-params.y(index))>min_real_delta_y)
        if ((current_max_y-current_min_y)>min_real_delta_y)
            % Yes - now, if this is the first extreme, check that you've
            % gone the requisite distance above the first point
            if ((last_extreme==0)& ...
                    ((current_max_y-params.y(1))>min_real_delta_y)) | ...
                (last_extreme==-1)
                % Are we far enough away from the maxium to know for sure?
                if ((last_extreme==0) & ...
                        (index>params.min_x_index_spacing)) | ...
                    ((last_extreme==-1) && ...
                        (index-current_max_index) > ...
                            params.min_x_index_spacing)
                    % We've defined a maximum
                    max_indices=[max_indices current_max_index];
                    current_min_y=params.y(index);
                    current_min_index=index;
                    current_max_y=params.y(index);
                    current_max_index=index;
                    last_extreme=1;
                end
            end
        end
    end
end
        
% Set return data
peak_data.min_indices=min_indices;
peak_data.max_indices=max_indices;    

            
    

