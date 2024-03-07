classdef ImFInfoDialog_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        GelImageFileInformationUIFigure  matlab.ui.Figure
        Tree  matlab.ui.container.Tree
    end

    
    properties (Access = private)
        GelBoxApp % Description
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, caller)
            
            movegui(app.GelImageFileInformationUIFigure,'center')

            app.GelBoxApp = caller;
            struct = caller.gel_data.image.imfinfo;
            parent = app.Tree;
            names = fieldnames(struct);

            for i = 1 : numel(names)
                if isa(struct.(names{i}),'char')
                    ch_node_text = sprintf('%s: %s', ...
                        convertCharsToStrings(names(i)), ...
                        convertCharsToStrings(struct.(names{i})));
                elseif isa(struct.(names{i}),'double')
                    [m,n] = size(struct.(names{i}));
                    if isempty(struct.(names{i}))
                        ch_node_text = sprintf('%s: [ ]', ...
                            convertCharsToStrings(names(i)));
                    elseif m == 1 && n==1
                        if rem(struct.(names{i}),1) == 0
                            formatting = '%s: %i';
                        else
                            formatting = '%s: %.3f';
                        end
                        ch_node_text = sprintf(formatting, ...
                            convertCharsToStrings(names(i)),struct.(names{i}));
                    elseif n > 4
                        ch_node_text = sprintf('%s: [%ix%i]', ...
                            convertCharsToStrings(names(i)),m,n);
                    else
                        formatting = '%s: [';
                        for j = 1:n
                            formatting = strcat(formatting,'%i, ');
                        end
                        formatting(end) = [];
                        formatting = strcat(formatting,']');
                        ch_node_text = sprintf(formatting, ...
                            convertCharsToStrings(names(i)), ...
                            struct.(names{i})(1:end));
                    end
                end

                uitreenode(parent,'Text',ch_node_text);
            end
            expand(app.Tree)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create GelImageFileInformationUIFigure and hide until all components are created
            app.GelImageFileInformationUIFigure = uifigure('Visible', 'off');
            app.GelImageFileInformationUIFigure.Position = [100 100 375 480];
            app.GelImageFileInformationUIFigure.Name = 'Gel Image File Information';

            % Create Tree
            app.Tree = uitree(app.GelImageFileInformationUIFigure);
            app.Tree.Position = [14 10 353 461];

            % Show the figure after all components are created
            app.GelImageFileInformationUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ImFInfoDialog_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.GelImageFileInformationUIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.GelImageFileInformationUIFigure)
        end
    end
end