classdef LayoutTableWindow_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        LaneLayoutUIFigure  matlab.ui.Figure
        UITable             matlab.ui.control.Table
    end


    properties (Access = private)
        GelBoxApp % Description
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, caller)
            movegui(app.LaneLayoutUIFigure,'center')
            app.GelBoxApp = caller;             
            app.UITable.Data = app.GelBoxApp.gel_data.layout.layout_table;
            app.UITable.ColumnName = app.GelBoxApp.gel_data.layout.layout_table.Properties.VariableNames;

        end

        % Close request function: LaneLayoutUIFigure
        function LaneLayoutUIFigureCloseRequest(app, event)
            delete(app)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create LaneLayoutUIFigure and hide until all components are created
            app.LaneLayoutUIFigure = uifigure('Visible', 'off');
            app.LaneLayoutUIFigure.Position = [100 100 630 421];
            app.LaneLayoutUIFigure.Name = 'Lane Layout';
            app.LaneLayoutUIFigure.CloseRequestFcn = createCallbackFcn(app, @LaneLayoutUIFigureCloseRequest, true);

            % Create UITable
            app.UITable = uitable(app.LaneLayoutUIFigure);
            app.UITable.ColumnName = {''; ''; ''; ''};
            app.UITable.RowName = {};
            app.UITable.Position = [16 11 603 394];

            % Show the figure after all components are created
            app.LaneLayoutUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = LayoutTableWindow_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.LaneLayoutUIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.LaneLayoutUIFigure)
        end
    end
end