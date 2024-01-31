classdef FittingOptionsDialog_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        FittingParametersUIFigure  matlab.ui.Figure
        UITable                    matlab.ui.control.Table
        UpdateFittingButton        matlab.ui.control.Button
    end


    properties (Access = private)
        GelBoxApp % Description
        n
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, caller)
            movegui(app.FittingParametersUIFigure,'center')
            app.GelBoxApp = caller;
            app.n = str2num(caller.BoxSelectionDropDown.Value);
           
            par_fields = {'Peak Location';'Amplitude';'Width';'Skew'};
            band_no_fields = {'1';'';'';''};

            no_of_bands = numel(caller.gel_data.fitting.par_est(app.n).band_no);
            
            mt.band_no_fields = repmat(band_no_fields,no_of_bands,1);
            m = 1;
            for i = 1:no_of_bands
                mt.band_no_fields{m} = sprintf('%i',i);
                m = m + numel(par_fields);
            end
            mt.par_fields = repmat(par_fields,no_of_bands,1);
            
            if no_of_bands > 1
                for i = 2:no_of_bands
                mt.par_fields{4*i-1,1} = 'Width Offset';
                mt.par_fields{4*i,1} = 'Skew Offset';
                end
            end

            names = fieldnames(caller.gel_data.fitting.par_est);
            par_names = names(2:end,:);

            l = 1;
            for u = 1 : no_of_bands
                for i = 1 : numel(par_names)
                    mt.par_est_values(l,1) = caller.gel_data.fitting.par_est(app.n).(par_names{i})(u);
                    mt.par_cal_values(l,1) = caller.gel_data.fitting.par_fit(app.n).(par_names{i})(u);
                    mt.par_con_values(l,1) = caller.gel_data.fitting.par_con(app.n).(par_names{i})(u);
                    l = l + 1;
                end
            end    
            
            t = struct2table(mt);

            app.UITable.Data = t;

        end

        % Button pushed function: UpdateFittingButton
        function UpdateFittingButtonPushed(app, event)
            UpdateFittingOptions(app.GelBoxApp);

            delete(app)
        end

        % Cell edit callback: UITable
        function StartingParametersTableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            new_param = table2struct(app.UITable.Data);
            names = fieldnames(app.GelBoxApp.gel_data.fitting.par_est);
            par_names = names(2:5,:);
            no_of_bands = numel(app.GelBoxApp.gel_data.fitting.par_est(app.n).band_no);            
            ll = 1;
            for ii = 1 : no_of_bands
                for jj = 1 : numel(par_names)  
                    app.GelBoxApp.gel_data.fitting.par_est(app.n).(par_names{jj})(ii) = new_param(ll,1).par_est_values;
                    app.GelBoxApp.gel_data.fitting.par_con(app.n).(par_names{jj})(ii) = new_param(ll,1).par_con_values;
                    ll = ll + 1;
                end                
            end
            app.GelBoxApp.gel_data.par_update(app.n) = 1;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create FittingParametersUIFigure and hide until all components are created
            app.FittingParametersUIFigure = uifigure('Visible', 'off');
            app.FittingParametersUIFigure.Position = [100 100 630 453];
            app.FittingParametersUIFigure.Name = 'Fitting Parameters';

            % Create UpdateFittingButton
            app.UpdateFittingButton = uibutton(app.FittingParametersUIFigure, 'push');
            app.UpdateFittingButton.ButtonPushedFcn = createCallbackFcn(app, @UpdateFittingButtonPushed, true);
            app.UpdateFittingButton.Position = [251 13 100 22];
            app.UpdateFittingButton.Text = 'Update Fitting';

            % Create UITable
            app.UITable = uitable(app.FittingParametersUIFigure);
            app.UITable.ColumnName = {'Band No'; 'Parameter'; 'Starting Parameter Estimate'; 'Calculated Parameter'; 'Constrain'};
            app.UITable.RowName = {};
            app.UITable.ColumnEditable = [false false true false true];
            app.UITable.CellEditCallback = createCallbackFcn(app, @StartingParametersTableCellEdit, true);
            app.UITable.Position = [14 45 599 393];

            % Show the figure after all components are created
            app.FittingParametersUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = FittingOptionsDialog_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.FittingParametersUIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.FittingParametersUIFigure)
        end
    end
end