classdef SelectedBoxDialog_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        SelectedBoxInformationUIFigure  matlab.ui.Figure
        Tree                            matlab.ui.container.Tree
    end

    
    properties (Access = private)
        GelBoxApp % Description
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, caller)
            
            movegui(app.SelectedBoxInformationUIFigure,'center')
            
            app.GelBoxApp = caller;
            struct = caller.gel_data.d_box.box;
            for i = 1 : length(struct)

            par_node_text = sprintf('Box %i',i);
            parent = uitreenode(app.Tree, 'Text',par_node_text);
            
            ch_node_text = sprintf('Box Position: [%.3f %.3f %.3f %.3f]',struct(i).position(1), ...
                struct(i).position(2),struct(i).position(3),struct(i).position(4));
            uitreenode(parent,'Text',ch_node_text);

            [m,n] = size(struct(i).inset);
            ch_node_text = sprintf('Box Size (Pixels): %d by %d',m,n);
            uitreenode(parent,'Text',ch_node_text);

            ch_node_text = sprintf('Total Area (Pixels): %.3f',struct(i).total_area);
            uitreenode(parent,'Text',ch_node_text);

            ch_node_text = sprintf('Background Area (Pixels): %.3f',struct(i).background_area);
            uitreenode(parent,'Text',ch_node_text);

            ch_node_text = sprintf('Number of Bands: %d',struct(i).fitting_mode);
            uitreenode(parent,'Text',ch_node_text);

            for j = 1:numel(struct(i).band_area)
                ch_node_text = sprintf('Band %i Area (Pixels): %.3f',j,struct(i).band_area(j));
                uitreenode(parent,'Text',ch_node_text);
            end
            expand(app.Tree)
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create SelectedBoxInformationUIFigure and hide until all components are created
            app.SelectedBoxInformationUIFigure = uifigure('Visible', 'off');
            app.SelectedBoxInformationUIFigure.Position = [100 100 375 480];
            app.SelectedBoxInformationUIFigure.Name = 'Selected Box Information';
            app.SelectedBoxInformationUIFigure.WindowStyle = 'modal';

            % Create Tree
            app.Tree = uitree(app.SelectedBoxInformationUIFigure);
            app.Tree.Position = [12 14 353 459];

            % Show the figure after all components are created
            app.SelectedBoxInformationUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SelectedBoxDialog_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.SelectedBoxInformationUIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.SelectedBoxInformationUIFigure)
        end
    end
end