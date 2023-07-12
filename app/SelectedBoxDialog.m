
classdef SelectedBoxDialog < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        SelectedBoxInformationUIFigure  matlab.ui.Figure
        Tree                            matlab.ui.container.Tree
        Node                            matlab.ui.container.TreeNode
        Node2uttut                      matlab.ui.container.TreeNode
        Node3                           matlab.ui.container.TreeNode
        Node4                           matlab.ui.container.TreeNode
        TextArea                        matlab.ui.control.TextArea
    end

    
    properties (Access = private)
        GelBoxApp % Description
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, caller)
            writelines(evalc('type(mfilename(''fullpath'')+".mlapp")'),mfilename('fullpath')+".m");

            % app.GelBoxApp = caller;
            % app.TextArea.Value = printstruct(app.GelBoxApp.gel_data.d_box, 'maxarray', 100);
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

            % Create TextArea
            app.TextArea = uitextarea(app.SelectedBoxInformationUIFigure);
            app.TextArea.Editable = 'off';
            app.TextArea.Position = [-439 -9 356 444];

            % Create Tree
            app.Tree = uitree(app.SelectedBoxInformationUIFigure);
            app.Tree.Position = [12 14 355 448];

            % Create Node
            app.Node = uitreenode(app.Tree);
            app.Node.Text = 'Node';

            % Create Node2uttut
            app.Node2uttut = uitreenode(app.Node);
            app.Node2uttut.Text = 'Node2: uttut';

            % Create Node3
            app.Node3 = uitreenode(app.Node);
            app.Node3.Text = 'Node3';

            % Create Node4
            app.Node4 = uitreenode(app.Node);
            app.Node4.Text = 'Node4';

            % Show the figure after all components are created
            app.SelectedBoxInformationUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SelectedBoxDialog(varargin)

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

