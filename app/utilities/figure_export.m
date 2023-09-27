function figure_export(varargin)
% This is essentially a wrapper for export_fig
% For more details on these options see
% http://sites.google.com/site/oliverwoodford/software/export_fig

params.output_file_string='';
params.output_type='eps';
params.renderer='painters';
params.dpi=250;
params.colorspace='rgb';
params.back_mode=[];
params.quality=100;

params=parse_pv_pairs(params,varargin);

% If output_type is empty do nothing
if (isempty(params.output_type))
    return;
end

% Generate output_file_string
ofs=[params.output_file_string '.' sprintf('%s',params.output_type)];

% Check whether the output file exists
if (exist(ofs))
    delete(ofs);
    pause(1);
end

% Code

% Work out where this function was called from
st=dbstack('-completenames');
calling_function=st(2).file;

% Check if we are sending out as svg
if (strcmp(params.output_type,'svg'))
    plot2svg(ofs,gcf);
    return
end    

% Do some clever checking for the the back_mode
if ((strcmp(params.output_type,'eps'))&(isempty(params.back_mode)))
    params.back_mode=1;
end
if ((strcmp(params.output_type,'pdf'))&(isempty(params.back_mode)))
    params.back_mode=1;
end
if ((strcmp(params.output_type,'png'))&(isempty(params.back_mode)))
    params.back_mode=1;
end
if ((strcmp(params.output_type,'tif'))&(isempty(params.back_mode)))
    params.back_mode=1;
end

dfc=get(gcf,'Color');
if (params.back_mode)
    dfc=get(gcf,'Color');
    set(gcf,'Color',[1 1 1]);
end
if ((strcmp(params.renderer,'zbuffer'))|| ...
        (strcmp(params.renderer,'opengl')))
    set(gcf,'Color','w');
end

if (~strcmp(params.output_type,'emf'))
    export_fig(sprintf('%s',ofs), ...
        sprintf('-r%d',params.dpi), ...
        sprintf('-%s',params.renderer), ...
        sprintf('-%s',params.colorspace), ...
        sprintf('-q%d',params.quality), ...
        '-nocrop');
else
    saveas(gcf,ofs);
end

if (strcmp(params.output_type,'eps'))
    disp('Correcting fonts');
    correct_fonts([params.output_file_string]);
end

if (params.back_mode)
    set(gcf,'Color',dfc);
end
if ((strcmp(params.renderer,'zbuffer'))|| ...
        (strcmp(params.renderer,'opengl')))
    set(gcf,'Color',dfc);
end

% Add in information about how this function was created
if (strcmp(params.output_type,'eps'))
    in_file=fopen(ofs,'r');
    out_file_string=[ofs 'temp'];
    out_file=fopen(out_file_string,'w');
    while (~feof(in_file))
        line_string=fgetl(in_file);
        if (strfind(line_string,'%%Title:'))
            fprintf(out_file,'%%%%Created by kens_export.m and %s\n', ...
                 calling_function);
        else
            fprintf(out_file,'%s\n',line_string);
        end
     end
     fclose(in_file);
     fclose(out_file);
     delete(ofs);
     copyfile(out_file_string,ofs,'f');
     delete(out_file_string);
end

% Png
if (strcmp(params.output_type,'png'))
    im=imread(ofs);
    imwrite(im,ofs, ...
        'Comment',sprintf('Created by kens_export.m and %s',calling_function));
end

% Tif
if (strcmp(params.output_type,'tif'))
    im=imread(ofs);
    imwrite(im,ofs, ...
        'Description',sprintf('Created by kens_export.m and %s',calling_function));
end
