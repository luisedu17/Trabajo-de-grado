classdef SC_Comparador_final < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                     matlab.ui.Figure
        GridLayout                   matlab.ui.container.GridLayout
        GridLayout10                 matlab.ui.container.GridLayout
        PanelBottom                  matlab.ui.container.Panel
        GridLayout8                  matlab.ui.container.GridLayout
        UITableSwitch                matlab.ui.control.Table
        UITableCapacitor             matlab.ui.control.Table
        PanelTop                     matlab.ui.container.Panel
        GridLayout3                  matlab.ui.container.GridLayout
        GridLayout9                  matlab.ui.container.GridLayout
        Lamp                         matlab.ui.control.Lamp
        BidirectionalSwitchedCapacitorConverterAnalyzerLabel  matlab.ui.control.Label
        GridLayout2                  matlab.ui.container.GridLayout
        GridLayout5                  matlab.ui.container.GridLayout
        TabGroup                     matlab.ui.container.TabGroup
        TimeResponseTab              matlab.ui.container.Tab
        UIAxes                       matlab.ui.control.UIAxes
        EfficiencyGraphTab           matlab.ui.container.Tab
        UIAxes2                      matlab.ui.control.UIAxes
        PerformanceParametersTab     matlab.ui.container.Tab
        GridLayout11                 matlab.ui.container.GridLayout
        Panel_8                      matlab.ui.container.Panel
        GridLayout13                 matlab.ui.container.GridLayout
        UITable1_3                   matlab.ui.control.Table
        Panel_7                      matlab.ui.container.Panel
        GridLayout12                 matlab.ui.container.GridLayout
        UITable3_1                   matlab.ui.control.Table
        GridLayout6                  matlab.ui.container.GridLayout
        TreeDickson                  matlab.ui.container.CheckBoxTree
        DicksonNode                  matlab.ui.container.TreeNode
        VoltajesNode                 matlab.ui.container.TreeNode
        DckVin                       matlab.ui.container.TreeNode
        DckVout                      matlab.ui.container.TreeNode
        CurrentsNode_3               matlab.ui.container.TreeNode
        IinNode_3                    matlab.ui.container.TreeNode
        IoutNode_3                   matlab.ui.container.TreeNode
        TreeLadder                   matlab.ui.container.CheckBoxTree
        LadderNode                   matlab.ui.container.TreeNode
        VoltagesNode_2               matlab.ui.container.TreeNode
        LadVin                       matlab.ui.container.TreeNode
        LadVout                      matlab.ui.container.TreeNode
        CurrentsNode_2               matlab.ui.container.TreeNode
        IinNode_2                    matlab.ui.container.TreeNode
        IoutNode_2                   matlab.ui.container.TreeNode
        TreeSP                       matlab.ui.container.CheckBoxTree
        SeriesParallelNode           matlab.ui.container.TreeNode
        VoltagesNode                 matlab.ui.container.TreeNode
        SPVin                        matlab.ui.container.TreeNode
        SPVout                       matlab.ui.container.TreeNode
        CurrentsNode                 matlab.ui.container.TreeNode
        IinNode                      matlab.ui.container.TreeNode
        IoutNode                     matlab.ui.container.TreeNode
        PanelLeft                    matlab.ui.container.Panel
        GridLayout4                  matlab.ui.container.GridLayout
        Panel_3                      matlab.ui.container.Panel
        MaximumPowerWEditField       matlab.ui.control.NumericEditField
        MaximumPowerWEditFieldLabel  matlab.ui.control.Label
        GridLayout7                  matlab.ui.container.GridLayout
        ClearButton                  matlab.ui.control.Button
        ButtonGroup                  matlab.ui.container.ButtonGroup
        Button1_3                    matlab.ui.control.RadioButton
        Button3_1                    matlab.ui.control.RadioButton
        ButtonRun                    matlab.ui.control.Button
        Panel_6                      matlab.ui.container.Panel
        SimulationTimesEditField     matlab.ui.control.NumericEditField
        SimulationTimesLabel         matlab.ui.control.Label
        Panel_5                      matlab.ui.container.Panel
        FrequencyKHzEditField        matlab.ui.control.NumericEditField
        FrequencyKHzEditFieldLabel   matlab.ui.control.Label
        Panel_4                      matlab.ui.container.Panel
        EfficiencyEditField          matlab.ui.control.NumericEditField
        EfficiencyEditFieldLabel     matlab.ui.control.Label
        Panel_2                      matlab.ui.container.Panel
        VoltageInputVEditField       matlab.ui.control.NumericEditField
        VoltageInputVEditFieldLabel  matlab.ui.control.Label
    end


    % Public properties that correspond to the Simulink model
    properties (Access = public, Transient)
        Simulation simulink.Simulation
    end

    
    properties (Access = private)
        Voltage_input double; %Input Voltage
        Maximum_Power_out double; %Maximum Power Output 
        Efficiency double; %Percentage Efficiency
        Frequency double; %Switching Frequency
        Sim_Time double; % Transient Timulation Time
        %Conversion_ratio double; %Conversion ratio
        Operation_Type string; %Conversion type to accsess struct
        Results struct = struct(); % Structure to store results
        Component_values struct = struct(); % Structure for component values
        Charge_Vector_C struct = struct(); % Charge vector for configuration C
        Charge_Vector_SW struct = struct(); % Charge vector for configuration SW
        Mat_Topologys struct = struct(); % Matrices for topologies
        Netlist_values struct = struct(); % Values for the netlist
        Ports struct = struct(); % Port connections
        %% Signals
        Temp_Results struct = struct();
        time struct = struct();
        SP_Vin struct = struct();
        SP_Vout struct = struct();
        SP_Iin struct = struct();
        SP_Iout struct = struct();
        Lad_Vin struct = struct();
        Lad_Vout struct = struct();
        Lad_Iin struct = struct();
        Lad_Iout struct = struct();
        Dck_Vin struct = struct();
        Dck_Vout struct = struct();
        Dck_Iin struct = struct();
        Dck_Iout struct = struct();
        ResultsAvg struct = struct();
        ResultsAvgtem struct = struct();
        Efficiency_Fsw struct = struct();
        Performance_Params struct = struct();
        

    end
    methods (Access = private)
        function Avg = Avg_value(app, Freq, time_vector, signal_vector)
            % Switching period
            T  = 1/Freq;                 
            % Number of cycles to average
            N  = 10;                    
            % Duration of the averaging window
            T_win = N*T;                
            % Final time of the simulation
            t_end = time_vector(end);   
            % Start time of the window (avoid negative values)
            t_start = max(time_vector(1), t_end - T_win);  
    
            % Logical index for samples inside the window
            idx = (time_vector >= t_start & time_vector <= t_end);
            % Time vector within the window
            t_seg = time_vector(idx);
            % Signal values within the window
            v_seg = signal_vector(idx);
    
            % Time-weighted average using trapezoidal integration
            Avg = trapz(t_seg, v_seg) / (t_seg(end) - t_seg(1));
        end
        function computeAverages(app, currentOp, fsw)
            %Valores promedio condisiones de diseño 
            topologies = {'SP', 'Lad', 'Dck'};
            signals    = {'Vout', 'Iout', 'Iin'};
    
            for t=1:length(topologies)
                topo=topologies{t};
                for s=1:length(signals)
                    sig=signals{s};
                    time_vector=app.time.(currentOp);
                    signal_vector=app.(sprintf('%s_%s',topo,sig)).(currentOp);
                    Avg = Avg_value(app, fsw, time_vector, signal_vector);
                    app.ResultsAvg.(sprintf('%s_%s',topo,sig)).(currentOp)=Avg;
                end
            end
        end
        function computeAveragestemp(app, currentOp, fsw)
            %Valores promedio condisiones para R_line R_load
            topologies = {'SP', 'Lad', 'Dck'};
            signals    = {'Vout', 'Iout', 'Iin'};
            levels = {'Vlow','Vhigh','RLhigh'};
    
            for t=1:length(topologies)
                topo=topologies{t};
                for s=1:length(signals)
                    sig=signals{s};
                    for l=1:length(levels)
                        lev=levels{l};
                        if isfield(app.Temp_Results.(sprintf('%s_%s',topo,sig)).(currentOp), lev)
                            time_vector = app.Temp_Results.time.(currentOp).(lev);
                            signal_vector  = app.Temp_Results.(sprintf('%s_%s',topo,sig)).(currentOp).(lev);
                            Avg = Avg_value(app, fsw, time_vector, signal_vector);
                            app.ResultsAvgtem.(sprintf('%s_%s',topo,sig)).(currentOp).(lev) = Avg;
                        end
                    end
                end
            end
        end
        function computeAveragesFreq(app, currentOp, fsw)
            %Valores promedio para grafivar eff vs freq
            Freq = [0.6*fsw, 0.8*fsw, fsw, 1.2*fsw, 1.4*fsw];
            topologies = {'SP','Lad','Dck'};
            signals    = {'Vout','Iout','Iin'};
            for f=1:length(Freq)
                fr =Freq(f);
                fieldf = sprintf('V_%.2f', fr); 
                fswField = strrep(fieldf,'.','_');
                for t = 1:length(topologies)
                    topo = topologies{t};
                    for s = 1:length(signals)
                        sig = signals{s};
                        if isfield(app.Temp_Results.(sprintf('%s_%s',topo,sig)).(currentOp), fswField)
                            time_vector = app.Temp_Results.time.(currentOp).(fswField);
                            signal_vector  = app.Temp_Results.(sprintf('%s_%s',topo,sig)).(currentOp).(fswField);
                            Avg = Avg_value(app, fr, time_vector, signal_vector);
                            app.ResultsAvgtem.(sprintf('%s_%s',topo,sig)).(currentOp).(fswField) = Avg;
                        end
                    end
                end
            end
        end

        function EfficiencyFreq (app, currentOp, fsw, Vin)
            Topologies = {'SP', 'Lad', 'Dck'};
            for t = 1:length(Topologies)
                topo = Topologies{t};
                if isfield(app.Efficiency_Fsw, topo) && isfield(app.Efficiency_Fsw.(topo), currentOp)
                    app.Efficiency_Fsw.(topo).(currentOp) = struct();
                end
            end
            Freq = [0.6*fsw, 0.8*fsw, fsw, 1.2*fsw, 1.4*fsw];
            for f=1:length(Freq)
                fr= Freq(f);
                fieldf = sprintf('V_%.2f', fr);  
                fswField = strrep(fieldf,'.','_');
                for t=1:length(Topologies)
                    topo = Topologies{t};
                    Vout_avg = app.ResultsAvgtem.(sprintf('%s_Vout',topo)).(currentOp).(fswField);
                    Iout_avg = app.ResultsAvgtem.(sprintf('%s_Iout',topo)).(currentOp).(fswField);
                    Iin_avg  = app.ResultsAvgtem.(sprintf('%s_Iin',topo)).(currentOp).(fswField);
                    Pout = Vout_avg * Iout_avg;
                    Pin  = Vin * Iin_avg;
                    effi_fsw  = Pout / Pin;
                    app.Efficiency_Fsw.(topo).(currentOp).(fswField) = abs(effi_fsw);
                end
            end
            
        end

        function PerformanceParams(app, Vin,currentOp)
            PerformanceP = {'Vout', 'R_line', 'R_load', 'P_loss', 'Eff'};
            Topologies = {'SP', 'Lad', 'Dck'};
            for t=1:length(Topologies)
                topo = Topologies{t};
                for p=1:length(PerformanceP)
                    perfor = PerformanceP{p};
                    switch perfor
                        case 'Vout'
                            app.Performance_Params.(perfor).(topo).(currentOp)=app.ResultsAvg.(sprintf('%s_Vout',topo)).(currentOp);
                        case 'R_line'
                            %'Vlow','Vhigh','RLhigh'
                            Vout2=app.ResultsAvgtem.(sprintf('%s_Vout',topo)).(currentOp).Vhigh;
                            Vout1=app.ResultsAvgtem.(sprintf('%s_Vout',topo)).(currentOp).Vlow;
                            Vin2=1.1*Vin;
                            Vin1=0.9*Vin;
                            Rline=((Vout2-Vout1)/(Vin2-Vin1));
                            app.Performance_Params.(perfor).(topo).(currentOp)=Rline;
                        case 'R_load'
                            Vouti2=app.ResultsAvg.(sprintf('%s_Vout',topo)).(currentOp);
                            Vouti1=app.ResultsAvgtem.(sprintf('%s_Vout',topo)).(currentOp).RLhigh;
                            I1=app.ResultsAvgtem.(sprintf('%s_Iout',topo)).(currentOp).RLhigh;
                            I2=app.ResultsAvg.(sprintf('%s_Iout',topo)).(currentOp);
                            Rload=((Vouti2-Vouti1)/(I2-I1));
                            app.Performance_Params.(perfor).(topo).(currentOp)=Rload;
                        case 'P_loss'
                            Vout=app.ResultsAvg.(sprintf('%s_Vout',topo)).(currentOp);
                            Iout=app.ResultsAvg.(sprintf('%s_Iout',topo)).(currentOp);
                            Iin=app.ResultsAvg.(sprintf('%s_Iin',topo)).(currentOp);
                            ploss=abs((Vin*Iin))-abs((Vout*Iout));
                            app.Performance_Params.(perfor).(topo).(currentOp)=ploss;
                        case 'Eff'
                            Voute=app.ResultsAvg.(sprintf('%s_Vout',topo)).(currentOp);
                            Ioute=app.ResultsAvg.(sprintf('%s_Iout',topo)).(currentOp);
                            Iine=app.ResultsAvg.(sprintf('%s_Iin',topo)).(currentOp);
                            effi=(Voute*Ioute)/(Vin*Iine);
                            app.Performance_Params.(perfor).(topo).(currentOp)=abs(effi);
                    end
                end
            end
        end
        function updateUITable3_1(app,currentOp)
            % Topologías en el mismo orden que las filas de la tabla
            Topologies = {'Lad','Dck','SP'};   % <-- ojo: RowName de tu tabla
            Params     = {'Vout','R_line','R_load','P_loss','Eff','Mssl','Mfsl'}; % Columnas
        
            % Prealocar matriz de celdas (3 filas x 7 columnas)
            data = cell(length(Topologies), length(Params));
        
            % Llenar con los resultados
            for t = 1:length(Topologies)
                topo = Topologies{t};
                for p = 1:length(Params)
                    perfor = Params{p};
                    if isfield(app.Performance_Params.(perfor).(topo), currentOp)
                        value = app.Performance_Params.(perfor).(topo).(currentOp);
                        % Redondear o formatear si deseas
                        data{t,p} = round(value,4);
                    else
                        data{t,p} = NaN; % si no existe el campo
                    end
                end
            end
        
            % Asignar a la tabla
            if currentOp == "ot3_1"
                app.UITable3_1.Data = data;
            elseif currentOp == "ot1_3"
                app.UITable1_3.Data = data;
            else
                warning('No se reconoce la operación: %s', currentOp);
            end
            
        end
    end

    methods (Access = private)
        function SimulationParam(app, ParamType, Vin, Rlmin, Freq, currentOp,netlist, ltspice_path)
                originalNetlist = netlist;  % Guardar netlist original con todos los placeholders
                switch ParamType
                    case "Vin"
                        values = [0.9*Vin, Vin, 1.1*Vin];
                        label = '{Vin}';
                    case "Rlmin"
                        values = [5*Rlmin];
                        label = '{RLmin}';
                    case "Freq"
                        values = [0.6*Freq, 0.8*Freq, Freq, 1.2*Freq, 1.4*Freq];
                        label = '{Freq}';
                    otherwise
                        error('Parametro no soportado');
                end
                for i = 1:length(values)
                    tempNetlist = originalNetlist;   % Copia limpia cada iteración
                    tempValue = values(i);
                    tempNetlist = strrep(tempNetlist, label, num2str(tempValue));

                    % Reemplazar los otros parámetros fijos
                    switch ParamType
                        case "Vin"
                            tempNetlist = strrep(tempNetlist, '{RLmin}', num2str(Rlmin));
                            tempNetlist = strrep(tempNetlist, '{Freq}', num2str(Freq));
                        case "Rlmin"
                            tempNetlist = strrep(tempNetlist, '{Vin}', num2str(Vin));
                            tempNetlist = strrep(tempNetlist, '{Freq}', num2str(Freq));
                        case "Freq"
                            tempNetlist = strrep(tempNetlist, '{Vin}', num2str(Vin));
                            tempNetlist = strrep(tempNetlist, '{RLmin}', num2str(Rlmin));
                    end

                    fileName = sprintf('Temp_%s_%s_%.2f.net', currentOp, ParamType, tempValue);
                    fid = fopen(fileName, 'w');
                    fprintf(fid, '%s', tempNetlist);
                    fclose(fid);
                    %% Ruta al ejecutable de LTspice
                    %ltspice_path = '"C:\Program Files\ADI\LTspice\LTspice.exe"';  % Ajusta la ruta según tu instalación
                    %% Ruta al archivo de netlist modificado
                    netlist_path = fileName;

                    %% Comando para ejecutar LTspice en modo batch
                    command = [ltspice_path ' -b ' netlist_path];
                    system(command);
                    filedataName = sprintf('Temp_%s_%s_%.2f.raw', currentOp, ParamType, tempValue);
                    data = LTspice2Matlab(filedataName);
                    switch ParamType
                        case "Vin"
                            if tempValue == Vin
                                %% Extraer tiempo y señales
                                app.time.(currentOp) = data.time_vect;
                                %%Serie-Paralelo signals
                                app.SP_Vin.(currentOp)= data.variable_mat(strcmp(data.variable_name_list, 'V(spvin)'), :);
                                app.SP_Vout.(currentOp)= data.variable_mat(strcmp(data.variable_name_list, 'V(spvout)'), :);
                                app.SP_Iin.(currentOp)= data.variable_mat(strcmp(data.variable_name_list, 'I(V1)'), :);
                                app.SP_Iout.(currentOp)= data.variable_mat(strcmp(data.variable_name_list, 'I(RLSP)'), :);

                                %%Ladder signals
                                app.Lad_Vin.(currentOp)= data.variable_mat(strcmp(data.variable_name_list, 'V(ladvin)'), :);
                                app.Lad_Vout.(currentOp)= data.variable_mat(strcmp(data.variable_name_list, 'V(ladvout)'), :);
                                app.Lad_Iin.(currentOp)= data.variable_mat(strcmp(data.variable_name_list, 'I(V3)'), :);
                                app.Lad_Iout.(currentOp)= data.variable_mat(strcmp(data.variable_name_list, 'I(RLLAD)'), :);
                                %%Dickson signals
                                app.Dck_Vin.(currentOp)= data.variable_mat(strcmp(data.variable_name_list, 'V(dckvin)'), :);
                                app.Dck_Vout.(currentOp)= data.variable_mat(strcmp(data.variable_name_list, 'V(dckvout)'), :);
                                app.Dck_Iin.(currentOp)= data.variable_mat(strcmp(data.variable_name_list, 'I(V2)'), :);
                                app.Dck_Iout.(currentOp)= data.variable_mat(strcmp(data.variable_name_list, 'I(RLDCK)'), :);

                            elseif tempValue == 0.9*Vin
                                app.Temp_Results.time.(currentOp).Vlow=data.time_vect;
                                app.Temp_Results.SP_Vout.(currentOp).Vlow= data.variable_mat(strcmp(data.variable_name_list, 'V(spvout)'), :);
                                app.Temp_Results.Lad_Vout.(currentOp).Vlow= data.variable_mat(strcmp(data.variable_name_list, 'V(ladvout)'), :);
                                app.Temp_Results.Dck_Vout.(currentOp).Vlow= data.variable_mat(strcmp(data.variable_name_list, 'V(dckvout)'), :);
                            else
                                app.Temp_Results.time.(currentOp).Vhigh=data.time_vect;
                                app.Temp_Results.SP_Vout.(currentOp).Vhigh= data.variable_mat(strcmp(data.variable_name_list, 'V(spvout)'), :);
                                app.Temp_Results.Lad_Vout.(currentOp).Vhigh= data.variable_mat(strcmp(data.variable_name_list, 'V(ladvout)'), :);
                                app.Temp_Results.Dck_Vout.(currentOp).Vhigh= data.variable_mat(strcmp(data.variable_name_list, 'V(dckvout)'), :);
                            end
                        case "Rlmin"
                            app.Temp_Results.time.(currentOp).RLhigh=data.time_vect;
                            app.Temp_Results.SP_Vout.(currentOp).RLhigh= data.variable_mat(strcmp(data.variable_name_list, 'V(spvout)'), :);
                            app.Temp_Results.Lad_Vout.(currentOp).RLhigh= data.variable_mat(strcmp(data.variable_name_list, 'V(ladvout)'), :);
                            app.Temp_Results.Dck_Vout.(currentOp).RLhigh= data.variable_mat(strcmp(data.variable_name_list, 'V(dckvout)'), :);
                            app.Temp_Results.SP_Iout.(currentOp).RLhigh= data.variable_mat(strcmp(data.variable_name_list, 'I(RLSP)'), :);
                            app.Temp_Results.Lad_Iout.(currentOp).RLhigh= data.variable_mat(strcmp(data.variable_name_list, 'I(RLLAD)'), :);
                            app.Temp_Results.Dck_Iout.(currentOp).RLhigh= data.variable_mat(strcmp(data.variable_name_list, 'I(RLDCK)'), :);
                        case "Freq"
                            fieldf = sprintf('V_%.2f', tempValue); 
                            fsw = strrep(fieldf,'.','_'); 
                            app.Temp_Results.time.(currentOp).(fsw)=data.time_vect;
                            app.Temp_Results.SP_Vout.(currentOp).(fsw)= data.variable_mat(strcmp(data.variable_name_list, 'V(spvout)'), :);
                            app.Temp_Results.Lad_Vout.(currentOp).(fsw)= data.variable_mat(strcmp(data.variable_name_list, 'V(ladvout)'), :);
                            app.Temp_Results.Dck_Vout.(currentOp).(fsw)= data.variable_mat(strcmp(data.variable_name_list, 'V(dckvout)'), :);
                            app.Temp_Results.SP_Iout.(currentOp).(fsw)= data.variable_mat(strcmp(data.variable_name_list, 'I(RLSP)'), :);
                            app.Temp_Results.Lad_Iout.(currentOp).(fsw)= data.variable_mat(strcmp(data.variable_name_list, 'I(RLLAD)'), :);
                            app.Temp_Results.Dck_Iout.(currentOp).(fsw)= data.variable_mat(strcmp(data.variable_name_list, 'I(RLDCK)'), :);
                            app.Temp_Results.SP_Iin.(currentOp).(fsw)= data.variable_mat(strcmp(data.variable_name_list, 'I(V1)'), :);
                            app.Temp_Results.Lad_Iin.(currentOp).(fsw)= data.variable_mat(strcmp(data.variable_name_list, 'I(V3)'), :);
                            app.Temp_Results.Dck_Iin.(currentOp).(fsw)= data.variable_mat(strcmp(data.variable_name_list, 'I(V2)'), :);

                    end
                end
           
        end
    end
    
    methods (Access = private)
        function ComponentValues(app, Conversion_ratio, Effi, Topology, Operation_Type, Vin, R_ssl, F_sw)
            for i=1:length(Topology)
                Act_Topology=Topology(i);
                switch Act_Topology
                    case "Serie_Paralelo"
                        topo = 'SP';
                    case "Ladder"
                        topo = 'Lad';
                    case "Dickson"
                        topo = 'Dck'
                end
                %%
                %Data Upload
                
                Bc=app.Mat_Topologys.(Act_Topology).(Operation_Type).V_cap;
                Vsw=app.Mat_Topologys.(Act_Topology).(Operation_Type).V_sw;
                bin=app.Mat_Topologys.(Act_Topology).(Operation_Type).bin;
                % calculos
                V_circuit=-(Bc^-1)*bin*Vin;
                V=[Vin; V_circuit];
                V_circuit(end)=[];
                V_cap=V_circuit;
                V_r=Vsw*V;
                % Save Data 
                % -Vector de tensiones en los capacitores
                % -Vetor de tensiones en los Switch en OFF
                app.Results.(Act_Topology).V_cap= V_cap;
                app.Results.(Act_Topology).V_r= V_r;
                %%
                
                %Data Upload
                Charge_Vector_Ci=app.Charge_Vector_C.(Act_Topology).(Operation_Type);
                Charge_Vector_SWi=app.Charge_Vector_SW.(Act_Topology).(Operation_Type);
                Vci_Work=app.Results.(Act_Topology).V_cap;
                Vswi_Work=app.Results.(Act_Topology).V_r;
            
                %Calculando valores de componentes
                if length(Charge_Vector_Ci)==length(Vci_Work)
                    n_cap= length(Charge_Vector_Ci);
                    acu= 0;
                    for j=1:n_cap
                        acu=acu+abs(Charge_Vector_Ci(j)*Vci_Work(j));
                    end
                    app.Performance_Params.Mssl.(topo).(Operation_Type)=(2*(Conversion_ratio*Vin*Effi)^2)/acu^2;
                    E_tot=(1/(2*R_ssl*F_sw))*(acu)^2;
                    acu1=0;
                    for j=1:n_cap
                    acu1=acu1+abs(Charge_Vector_Ci(j)*Vci_Work(j));
                    end
                    Capacitors_value= zeros(1,n_cap);
                    for j=1:n_cap
                        Capacitors_value(j)=(abs(Charge_Vector_Ci(j))/abs(Vci_Work(j)))*((2*E_tot)/acu1);
                    end
                else
                    sprintf("The arrays are not the same length")
                end
            
                if length(Charge_Vector_SWi)==length(Vswi_Work)
                    R_fsl=R_ssl;
                    n_sw= length(Charge_Vector_SWi);
                    acu=0;
                    for j=1:n_sw
                        acu=acu+abs(Charge_Vector_SWi(j)*Vswi_Work(j));
                    end
                    app.Performance_Params.Mfsl.(topo).(Operation_Type)=((Conversion_ratio*Vin*Effi)^2)/(2*(acu^2));
                    Xtot=(2/R_fsl)*(acu)^2;
                    Conductance=zeros(1,n_sw);
                    Switch_Ron=zeros(1,n_sw);
                    for j=1:n_sw
                        Conductance(j)=abs(Charge_Vector_SWi(j)/Vswi_Work(j))*(Xtot/acu);
                        Switch_Ron(j)=1/Conductance(j);
                    end
                else
                    sprintf("The arrays are not the same length")
                end
                % Save Data 
                % -Vector de valores de capacitancia
                % -Vetor de Resistencia de encendido de los switchs
                app.Component_values.(Act_Topology).(Operation_Type).Capacitors= Capacitors_value;
                app.Component_values.(Act_Topology).(Operation_Type).Ron_switch= Switch_Ron;
                
            end
        end
    end


    methods (Access = private)
        
        function Graphics(app)

             hold(app.UIAxes, 'on'); % Mantener señales previas
             % Obtener nodos seleccionados de cada topología de forma segura
            selectedSP = app.getCheckedNodesText(app.TreeSP);
            selectedLad = app.getCheckedNodesText(app.TreeLadder);
            selectedDck = app.getCheckedNodesText(app.TreeDickson);
            % Mapeo de señales basado en la topología
            signalsMap = struct(...
                'Serie_Paralelo', struct('Vin', app.SP_Vin, 'Vout', app.SP_Vout, ...
                'Iin', app.SP_Iin, 'Iout', app.SP_Iout), ...
                'Ladder', struct('Vin', app.Lad_Vin, 'Vout', app.Lad_Vout, ...
                'Iin', app.Lad_Iin, 'Iout', app.Lad_Iout), ...
                'Dickson', struct('Vin', app.Dck_Vin, 'Vout', app.Dck_Vout, ...
                'Iin', app.Dck_Iin, 'Iout', app.Dck_Iout) ...
            );

            % Borrar solo las señales que ya no están seleccionadas
            delete(findall(app.UIAxes, 'Tag', 'DynamicPlot'));
            % Separar señales de voltaje y corriente
            voltSignals = {'Vin', 'Vout'}; % Definir cuáles son voltajes
            currentSignals = {'Iin', 'Iout'}; % Definir cuáles son corrientes

            colorMap = struct( ...
                'Serie_Paralelo', 'b', ... % Azul para SP
                'Ladder', 'r', ...         % Rojo para Ladder
                'Dickson', 'g' ...         % Verde para Dickson
            );
            fullNames  = struct('Serie_Paralelo','Series-Parallel', ...
                            'Ladder','Ladder', ...
                            'Dickson','Dickson');

            % Función anidada para evitar repetir código
            function plotSignals(tipo, selectedNodes)
                if isempty(selectedNodes)
                    return; % No hacer nada si no hay selección
                end
                for i = 1:numel(selectedNodes)
                    nodeText = selectedNodes{i}; % Nombre de la señal seleccionada
                    if isfield(signalsMap.(tipo), nodeText) % Verificar si existe
                        if any(strcmp(nodeText, voltSignals)) % Si es voltaje
                            yyaxis(app.UIAxes, 'left');
                            ylabel(app.UIAxes, 'Voltaje (V)');
                        elseif any(strcmp(nodeText, currentSignals)) % Si es corriente
                            yyaxis(app.UIAxes, 'right');
                            ylabel(app.UIAxes, 'Corriente (A)');
                        end
                        %lineColor = colorMap.(tipo);
                        plot(app.UIAxes, app.time.(app.Operation_Type), signalsMap.(tipo).(nodeText).(app.Operation_Type), ...
                             'DisplayName', sprintf('%s - %s', fullNames.(tipo), nodeText), 'Tag', 'DynamicPlot', ...
                             'LineStyle', '-','Color', colorMap.(tipo), 'Marker', 'none');
                    end
                end
            end

            % Graficar solo las señales seleccionadas
            plotSignals('Serie_Paralelo', selectedSP);
            plotSignals('Ladder', selectedLad);
            plotSignals('Dickson', selectedDck);
            legend(app.UIAxes, 'show', 'Location', 'best'); % Mostrar leyenda
            hold(app.UIAxes, 'off'); % Liberar hold
            
        end
        
        function selectedNodes = getCheckedNodesText(app, tree)
            if isempty(tree.CheckedNodes)
                selectedNodes = {}; % Si no hay nodos seleccionados, devuelve un array vacío
            else
                selectedNodes = {tree.CheckedNodes.Text}; % Extrae los textos de los nodos
            end
        end

        function plotEfficiencyFreq(app)
        % Limpia el eje
        freqs=[];
        effs=[];
        cla(app.UIAxes2);
        hold(app.UIAxes2, 'on');
        delete(findall(app.UIAxes2, 'Tag', 'EffPlot'));
    
        % Topologías y colores
        topologies = {'SP','Lad','Dck'};
        colorMap   = struct('SP','b','Lad','r','Dck','g');
        fullNames  = struct('SP','Series-Parallel', ...
                            'Lad','Ladder', ...
                            'Dck','Dickson');
    
        % Tipo de operación actual (3:1 o 1:3)
        currentOp = app.Operation_Type;
    
        for k = 1:length(topologies)
            topo = topologies{k};
    
            % Verificar que existan datos para esta topología y operación
            if isfield(app.Efficiency_Fsw, topo) && ...
               isfield(app.Efficiency_Fsw.(topo), currentOp)
    
                freqFields = fieldnames(app.Efficiency_Fsw.(topo).(currentOp));
                freqs = zeros(size(freqFields));
                effs  = zeros(size(freqFields));
    
                for i = 1:length(freqFields)
                    % Convertir el nombre del campo en valor de frecuencia
                    rawField   = freqFields{i};
                    numericStr = strrep(rawField(3:end), '_', '.');
                    freqs(i)   = str2double(numericStr);
                    effs(i)    = app.Efficiency_Fsw.(topo).(currentOp).(rawField);
                end
    
                % Ordenar por frecuencia para una curva ordenada
                [freqs, idx] = sort(freqs);
                effs = effs(idx);
    
                % Graficar con nombre completo en la leyenda
                plot(app.UIAxes2, freqs, effs, '-o', ...
                    'Color', colorMap.(topo), ...
                    'DisplayName', fullNames.(topo), ...
                    'Tag', 'EffPlot');
            end
        end
        legend(app.UIAxes2, 'show', 'Location', 'best');
        grid(app.UIAxes2, 'on');
        hold(app.UIAxes2, 'off');
    end
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.UIFigure.WindowState = "maximized";
            app.UITableCapacitor.RowName = {'Ladder','Dickson','Series-Parallel'};
            app.UITableCapacitor.Data=["", "","", "","";
                                       "", "","", "","";
                                       "", "","", "",""];
            app.UITableSwitch.RowName = {'Ladder','Dickson','Series-Parallel'};
            app.UITableSwitch.Data=["", "", "", "", "", "", "";
                                    "", "", "", "", "", "", "";
                                    "", "", "", "", "", "", ""];
            app.UITable3_1.RowName = {'Ladder','Dickson','Series-Parallel'};
            app.UITable3_1.Data=["", "", "", "", "","","";
                                 "", "", "", "", "","","";
                                 "", "", "", "", "","",""];
            app.UITable1_3.RowName = {'Ladder','Dickson','Series-Parallel'};
            app.UITable1_3.Data=["", "", "", "", "","","";
                                 "", "", "", "", "","","";
                                 "", "", "", "", "","",""];
            %% Axes 
            cla(app.UIAxes);
            cla(app.UIAxes2);
            %Left Y axis
            yyaxis(app.UIAxes,"left");
            app.UIAxes.YLabel.String = "Voltage (V)";
            app.UIAxes.YColor = [0 0.4 0.8]; % blue
            %Right Y axes
            yyaxis(app.UIAxes,"right");
            app.UIAxes.YLabel.String = "Ampers (A)";
            app.UIAxes.YColor = [0.8 0.2 0]; % Orange
            %Common axes 
            app.UIAxes.XLabel.String = 'Time (s)';
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            % Tittle axes
            title(app.UIAxes, 'Time Response - SC Converter');
            % Config axes
            app.UIAxes.FontName = 'Arial';
            app.UIAxes.FontSize = 12;
            app.UIAxes.Box = 'on';
            app.UIAxes.Toolbar.Visible = 'on';
            %indicador
            app.Lamp.Color = [1 0 0]; % Rojo
            app.Lamp.Visible = 'on';
            % limits axes
            ylim(app.UIAxes, 'auto');
            xlim(app.UIAxes, 'auto');
            % Establecer valor inicial
            app.ButtonGroup.SelectedObject = app.Button3_1; % Selecciona 3:1 por defecto
            app.Operation_Type = 'ot3_1'; % Valor inicial explícito
            %app.Conversion_ratio=1/3;
            app.Results = struct();
            app.Component_values = struct();
            %
            app.Charge_Vector_C.Serie_Paralelo.ot3_1=[1/3 1/3];
            app.Charge_Vector_C.Serie_Paralelo.ot1_3=[-1 -1];
            app.Charge_Vector_C.Ladder.ot3_1=[2/3 -1/3 1/3];
            app.Charge_Vector_C.Ladder.ot1_3=[-2 1 -1];
            app.Charge_Vector_C.Dickson.ot3_1=[-1/3 1/3];
            app.Charge_Vector_C.Dickson.ot1_3=[1 -1];
            %
            app.Charge_Vector_SW.Serie_Paralelo.ot3_1=[1/3 1/3 1/3 1/3 1/3 2/3 1/3];
            app.Charge_Vector_SW.Serie_Paralelo.ot1_3=[1 1 1 1 1 2 1];
            app.Charge_Vector_SW.Ladder.ot3_1=[1/3 1/3 1/3 1/3 2/3 2/3];
            app.Charge_Vector_SW.Ladder.ot1_3=[1 1 1 1 2 2];
            app.Charge_Vector_SW.Dickson.ot3_1=[1/3 1/3,1/3 1/3 1/3 1/3 1/3];
            app.Charge_Vector_SW.Dickson.ot1_3=[1 1 1 1 1 1 1];
            %
            app.Mat_Topologys.Serie_Paralelo.ot3_1.V_cap=[1 1 1;
                                                         -1 0 1;
                                                         -1 1 0];
            app.Mat_Topologys.Serie_Paralelo.ot3_1.V_sw=[1 0 0 -1;
                                                         1 -1 0 0;
                                                         1 0 -1 -1;
                                                         0 0 0 1;
                                                         1 -1 -1 0;
                                                         1 -1 0 -1;
                                                         0 1 -1 -1];
            app.Mat_Topologys.Serie_Paralelo.ot3_1.bin=[-1; 0; 0];
            app.Mat_Topologys.Serie_Paralelo.ot1_3.V_cap=[-1 -1 1;
                                                       0 1 0;
                                                      -1 1 0];
            app.Mat_Topologys.Serie_Paralelo.ot1_3.V_sw=[0 -1 0 1;
                                                     0 -1 0 1;
                                                     -1 0 -1 1;
                                                     0 0 1 0;
                                                     0 -1 -1 1;
                                                     -1 -1 0 1;
                                                     1 0 0 0];
            app.Mat_Topologys.Serie_Paralelo.ot1_3.bin=[-1; -1; 0];
            app.Mat_Topologys.Ladder.ot3_1.V_cap=[1 0 1 1;
                                              1 -1 0 0;
                                             -1 0 0 1;
                                              0 1 -1 0];
            app.Mat_Topologys.Ladder.ot3_1.V_sw=[1 -1 0 -1 0;
                                             1 0 -1 0 -1;
                                             0 1 -1 0 1;
                                             1 0 0 -1 -1;
                                             0 0 0 0 1;
                                             0 0 0 0 1];
            app.Mat_Topologys.Ladder.ot3_1.bin=[-1; 0; 0; 0];
            app.Mat_Topologys.Ladder.ot1_3.V_cap=[-1 0 -1 1;
                                              -1 1 0 0;
                                              1 0 0 0;
                                              0 -1 1 0];
            app.Mat_Topologys.Ladder.ot1_3.V_sw=[-1 0 -1 0 1;
                                             -1 0 -1 0 1;
                                              0 0 0 1 0;
                                              0 1 0 0 0;
                                              1 0 0 0 0;
                                              1 0 0 0 0];
            app.Mat_Topologys.Ladder.ot1_3.bin=[-1; 0; -1; 0];
            app.Mat_Topologys.Dickson.ot3_1.V_cap=[0 1 1;
                                              -1 0 1;
                                               1 -1 1];
            app.Mat_Topologys.Dickson.ot3_1.V_sw=[1 -1 0 -1;
                                              1 -1 0 0;
                                              0 0 1 -1;
                                              0 0 0 1;
                                              0 0 0 1;
                                              0 0 0 1;
                                              0 0 0 1];
            app.Mat_Topologys.Dickson.ot3_1.bin=[-1; 0; 0];
            app.Mat_Topologys.Dickson.ot1_3.V_cap=[0 -1 1;
                                               1 0 0;
                                              -1 1 0];
            app.Mat_Topologys.Dickson.ot1_3.V_sw=[1 1 0 -1;
                                              1 0 0 -1;
                                              1 0 -1 0;
                                              1 0 0 0;
                                              1 0 0 0;
                                              1 0 0 0;
                                              1 0 0 0];
            app.Mat_Topologys.Dickson.ot1_3.bin=[-1; -1; -1];
            %
            app.Netlist_values.Serie_Paralelo.Capacitors.C1="{C1ValSP}";
            app.Netlist_values.Serie_Paralelo.Capacitors.C2="{C2ValSP}";
            app.Netlist_values.Serie_Paralelo.SW_Types.SWT1="{RonSPSwT1}";
            app.Netlist_values.Serie_Paralelo.SW_Types.SWT2="{RonSPSwT2}";
            app.Netlist_values.Serie_Paralelo.SW_Types.SWT3="{RonSPSwT3}";
            app.Netlist_values.Ladder.Capacitors.C1="{C1ValLAD}";
            app.Netlist_values.Ladder.Capacitors.C2="{C2ValLAD}";
            app.Netlist_values.Ladder.Capacitors.C3="{C3ValLAD}";
            app.Netlist_values.Ladder.SW_Types.SWT1="{RonLADSwT1}";
            app.Netlist_values.Ladder.SW_Types.SWT2="{RonLADSwT2}";
            app.Netlist_values.Dickson.Capacitors.C1="{C1ValDCK}";
            app.Netlist_values.Dickson.Capacitors.C2="{C2ValDCK}";
            app.Netlist_values.Dickson.SW_Types.SWT1="{RonDCKSwT1}";
            app.Netlist_values.Dickson.SW_Types.SWT2="{RonDCKSwT2}";
            %
            app.Ports.Serie_Paralelo.in1="{SPIn1}";
            app.Ports.Serie_Paralelo.in2="{SPIn2}";
            app.Ports.Serie_Paralelo.Cout="{SPCout}";
            app.Ports.Ladder.in1="{LADIn1}";
            app.Ports.Ladder.in2="{LADIn2}";
            app.Ports.Ladder.Cout="{LADCout}";
            app.Ports.Dickson.in1="{DCKIn1}";
            app.Ports.Dickson.in2="{DCKIn2}";
            app.Ports.Dickson.Cout="{DCKCout}";


        end

        % Value changed function: VoltageInputVEditField
        function VoltageInputVEditFieldValueChanged(app, event)
            app.Voltage_input = app.VoltageInputVEditField.Value;
        end

        % Value changed function: MaximumPowerWEditField
        function MaximumPowerWEditFieldValueChanged(app, event)
            app.Maximum_Power_out = app.MaximumPowerWEditField.Value;
            
        end

        % Value changed function: EfficiencyEditField
        function EfficiencyEditFieldValueChanged(app, event)
            app.Efficiency = app.EfficiencyEditField.Value;
            
        end

        % Value changed function: FrequencyKHzEditField
        function FrequencyKHzEditFieldValueChanged(app, event)
            app.Frequency = app.FrequencyKHzEditField.Value;
            
        end

        % Value changed function: SimulationTimesEditField
        function SimulationTimesEditFieldValueChanged(app, event)
            app.Sim_Time = app.SimulationTimesEditField.Value*1e-6;
            
        end

        % Selection changed function: ButtonGroup
        function ButtonGroupSelectionChanged(app, event)
            selectedButton = app.ButtonGroup.SelectedObject;
            Convertion_type=selectedButton.Text;
            switch Convertion_type
                case '3:1'
                    app.Operation_Type="ot3_1";
                case '1:3'
                    app.Operation_Type="ot1_3";
            end

            plotEfficiencyFreq(app)
        end

        % Button pushed function: ButtonRun
        function ButtonRunPushed(app, event)
            app.Lamp.Color = [1 0 0]; % Rojo
            config = jsondecode(fileread('config.json'));
            ltspice_path = ['"', config.ltspice_path, '"']; % entre comillas por si hay espacios
            namelist = config.netlist_base;
            %Access to Topologys
            E_p = app.Efficiency;
            Effi=E_p/100;
            F_sw= app.Frequency*1000;%Frequency Switching in Hz
            Topology=["Serie_Paralelo", "Ladder", "Dickson"];
            Op_Type=["ot3_1", "ot1_3"];
            CoutVL=0;
            CoutVH=0;
            for i=1:length(Op_Type)
                currentOp= Op_Type(i);
                switch currentOp
                    case "ot3_1"
                        Conversion_ratio=1/3;
                        Vin=app.Voltage_input;
                        RL_min= (((Effi*Vin*Conversion_ratio)^2)/app.Maximum_Power_out);
                        Imax=(app.Maximum_Power_out)/(Conversion_ratio*Vin*Effi);
                        CoutVL=Imax/(F_sw*(Conversion_ratio*Vin*0.01));
                    case "ot1_3"
                        Conversion_ratio=3;
                        Vin=app.Voltage_input/3;
                        RL_min= (((Effi*Vin*Conversion_ratio)^2)/app.Maximum_Power_out);
                        Imax=(app.Maximum_Power_out)/(Conversion_ratio*Vin*Effi);
                        CoutVH=Imax/(F_sw*(Conversion_ratio*Vin*0.01));
                end
                Ro_max= RL_min*((1-Effi)/Effi);
                R_ssl=Ro_max/sqrt(2);
                %% Component Values 
                ComponentValues(app, Conversion_ratio, Effi, Topology, currentOp, app.Voltage_input, R_ssl, F_sw);
                % Table Values
                DataCapacitors = [ ...
                       [app.Component_values.Ladder.ot3_1.Capacitors, CoutVL, CoutVH]; ...
                       [app.Component_values.Dickson.ot3_1.Capacitors, "", "", ""]; ...
                       [app.Component_values.Serie_Paralelo.ot3_1.Capacitors, "", "", ""] ...
                       ];
                app.UITableCapacitor.Data= DataCapacitors;
                DataRon_switch = [ ...
                       [app.Component_values.Ladder.ot3_1.Ron_switch], ""; ...
                       [app.Component_values.Dickson.ot3_1.Ron_switch]; ...
                       [app.Component_values.Serie_Paralelo.ot3_1.Ron_switch] ...
                       ];
                app.UITableSwitch.Data= DataRon_switch;
                %% Read Netlist File
                %netlist = fileread('SC_Converter1.net');
                netlist = fileread(namelist);
                for k=1:length(Topology)
                    Act_Topology=Topology(k);
                    capacitores = fieldnames(app.Netlist_values.(Act_Topology).Capacitors);
                    for j=1:length(capacitores)
                        capacitor=capacitores{j};
                        netlist = strrep(netlist,app.Netlist_values.(Act_Topology).Capacitors.(capacitor) , ...
                        num2str(app.Component_values.(Act_Topology).(currentOp).Capacitors(j)));
                    end
                    Switch_Types = fieldnames(app.Netlist_values.(Act_Topology).SW_Types);
                    for j=1:length(Switch_Types)
                        TypeSW=Switch_Types{j};
                        Ron_SWType=unique(app.Component_values.(Act_Topology).(currentOp).Ron_switch, 'stable');
                        netlist = strrep(netlist,app.Netlist_values.(Act_Topology).SW_Types.(TypeSW) , ...
                        num2str(Ron_SWType(j)));
                    end
                end

                %Netlist Load Time
                netlist = strrep(netlist, '{Time}', num2str(app.Sim_Time));
                switch currentOp
                    case "ot3_1"
                        netlist =strrep(netlist,"{SPIn1}","SPVin");
                        netlist =strrep(netlist,"{SPIn2}","SPVout");
                        netlist =strrep(netlist,"{SPCout}", ...
                            num2str(CoutVL));
                        netlist =strrep(netlist,"{DCKIn1}","DCKVin");
                        netlist =strrep(netlist,"{DCKIn2}","DCKVout");
                        netlist =strrep(netlist,"{DCKCout}", ...
                            num2str(CoutVL));
                        netlist =strrep(netlist,"{LADIn1}","LADVin");
                        netlist =strrep(netlist,"{LADIn2}","LADVout");
                        netlist =strrep(netlist,"{LADCout}", ...
                            num2str(CoutVL));
                    case "ot1_3"
                        netlist =strrep(netlist,"{SPIn1}","SPVout");
                        netlist =strrep(netlist,"{SPIn2}","SPVin");
                        netlist =strrep(netlist,"{SPCout}", ...
                            num2str(CoutVH));
                        netlist =strrep(netlist,"{DCKIn1}","DCKVout");
                        netlist =strrep(netlist,"{DCKIn2}","DCKVin");
                        netlist =strrep(netlist,"{DCKCout}", ...
                            num2str(CoutVH));
                        netlist =strrep(netlist,"{LADIn1}","LADVout");
                        netlist =strrep(netlist,"{LADIn2}","LADVin");
                        netlist =strrep(netlist,"{LADCout}", ...
                            num2str(CoutVH));
                end
                ParamType=["Vin", "Rlmin", "Freq"];
                for j=1:length(ParamType)
                    Act_ParamType=ParamType(j);
                    SimulationParam(app, Act_ParamType, Vin, RL_min, F_sw, currentOp, netlist, ltspice_path);

                end
                computeAverages(app, currentOp, F_sw);
                computeAveragestemp(app, currentOp, F_sw);
                computeAveragesFreq(app, currentOp, F_sw);
                EfficiencyFreq (app, currentOp, F_sw, Vin);
                PerformanceParams(app, Vin,currentOp);
                updateUITable3_1(app,currentOp);
            end
            app.Lamp.Color = [0 1 0]; % Verde
        end

        % Callback function: TreeSP
        function TreeSPCheckedNodesChanged(app, event)
            Graphics(app); 
        end

        % Callback function: TreeLadder
        function TreeLadderCheckedNodesChanged(app, event)
            Graphics(app);
        end

        % Callback function: TreeDickson
        function TreeDicksonCheckedNodesChanged(app, event)
            Graphics(app);
        end

        % Button pushed function: ClearButton
        function ClearButtonPushed(app, event)
            startupFcn(app);
            app.Lamp.Color = [1 0 0]; % Rojo
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1x'};
            app.GridLayout.RowHeight = {'0.15x', '0.65x', '0.2x'};
            app.GridLayout.ColumnSpacing = 1;
            app.GridLayout.RowSpacing = 1;
            app.GridLayout.Padding = [2 2 2 2];

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.GridLayout);
            app.GridLayout2.ColumnWidth = {'0.3x', '0.7x'};
            app.GridLayout2.RowHeight = {'1x'};
            app.GridLayout2.ColumnSpacing = 0;
            app.GridLayout2.RowSpacing = 0;
            app.GridLayout2.Padding = [0 0 0 0];
            app.GridLayout2.Layout.Row = 2;
            app.GridLayout2.Layout.Column = 1;

            % Create PanelLeft
            app.PanelLeft = uipanel(app.GridLayout2);
            app.PanelLeft.TitlePosition = 'centertop';
            app.PanelLeft.Title = 'Design parameters';
            app.PanelLeft.Layout.Row = 1;
            app.PanelLeft.Layout.Column = 1;
            app.PanelLeft.FontName = 'Arial';
            app.PanelLeft.FontWeight = 'bold';

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.PanelLeft);
            app.GridLayout4.ColumnWidth = {'1x'};
            app.GridLayout4.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout4.RowSpacing = 2;
            app.GridLayout4.Padding = [2 2 2 2];

            % Create Panel_2
            app.Panel_2 = uipanel(app.GridLayout4);
            app.Panel_2.TitlePosition = 'centertop';
            app.Panel_2.Layout.Row = 1;
            app.Panel_2.Layout.Column = 1;

            % Create VoltageInputVEditFieldLabel
            app.VoltageInputVEditFieldLabel = uilabel(app.Panel_2);
            app.VoltageInputVEditFieldLabel.FontName = 'Arial';
            app.VoltageInputVEditFieldLabel.FontWeight = 'bold';
            app.VoltageInputVEditFieldLabel.Position = [8 11 107 22];
            app.VoltageInputVEditFieldLabel.Text = 'Voltage Input(V)';

            % Create VoltageInputVEditField
            app.VoltageInputVEditField = uieditfield(app.Panel_2, 'numeric');
            app.VoltageInputVEditField.Limits = [0 Inf];
            app.VoltageInputVEditField.ValueChangedFcn = createCallbackFcn(app, @VoltageInputVEditFieldValueChanged, true);
            app.VoltageInputVEditField.Position = [126 11 46 22];

            % Create Panel_4
            app.Panel_4 = uipanel(app.GridLayout4);
            app.Panel_4.TitlePosition = 'centertop';
            app.Panel_4.Layout.Row = 3;
            app.Panel_4.Layout.Column = 1;

            % Create EfficiencyEditFieldLabel
            app.EfficiencyEditFieldLabel = uilabel(app.Panel_4);
            app.EfficiencyEditFieldLabel.FontName = 'Arial';
            app.EfficiencyEditFieldLabel.FontWeight = 'bold';
            app.EfficiencyEditFieldLabel.Position = [8 11 108 22];
            app.EfficiencyEditFieldLabel.Text = 'Efficiency(%)';

            % Create EfficiencyEditField
            app.EfficiencyEditField = uieditfield(app.Panel_4, 'numeric');
            app.EfficiencyEditField.Limits = [0 100];
            app.EfficiencyEditField.ValueChangedFcn = createCallbackFcn(app, @EfficiencyEditFieldValueChanged, true);
            app.EfficiencyEditField.Position = [126 11 47 22];

            % Create Panel_5
            app.Panel_5 = uipanel(app.GridLayout4);
            app.Panel_5.TitlePosition = 'centertop';
            app.Panel_5.Layout.Row = 4;
            app.Panel_5.Layout.Column = 1;

            % Create FrequencyKHzEditFieldLabel
            app.FrequencyKHzEditFieldLabel = uilabel(app.Panel_5);
            app.FrequencyKHzEditFieldLabel.FontName = 'Arial';
            app.FrequencyKHzEditFieldLabel.FontWeight = 'bold';
            app.FrequencyKHzEditFieldLabel.Position = [8 12 107 22];
            app.FrequencyKHzEditFieldLabel.Text = 'Frequency(KHz)';

            % Create FrequencyKHzEditField
            app.FrequencyKHzEditField = uieditfield(app.Panel_5, 'numeric');
            app.FrequencyKHzEditField.Limits = [0 Inf];
            app.FrequencyKHzEditField.ValueDisplayFormat = '%5.3g';
            app.FrequencyKHzEditField.ValueChangedFcn = createCallbackFcn(app, @FrequencyKHzEditFieldValueChanged, true);
            app.FrequencyKHzEditField.Position = [126 11 46 22];

            % Create Panel_6
            app.Panel_6 = uipanel(app.GridLayout4);
            app.Panel_6.TitlePosition = 'centertop';
            app.Panel_6.Layout.Row = 5;
            app.Panel_6.Layout.Column = 1;

            % Create SimulationTimesLabel
            app.SimulationTimesLabel = uilabel(app.Panel_6);
            app.SimulationTimesLabel.FontName = 'Arial';
            app.SimulationTimesLabel.FontWeight = 'bold';
            app.SimulationTimesLabel.Position = [8 12 119 22];
            app.SimulationTimesLabel.Text = 'Simulation Time(µs)';

            % Create SimulationTimesEditField
            app.SimulationTimesEditField = uieditfield(app.Panel_6, 'numeric');
            app.SimulationTimesEditField.Limits = [0 Inf];
            app.SimulationTimesEditField.ValueDisplayFormat = '%5.3g';
            app.SimulationTimesEditField.ValueChangedFcn = createCallbackFcn(app, @SimulationTimesEditFieldValueChanged, true);
            app.SimulationTimesEditField.Position = [126 11 46 22];

            % Create GridLayout7
            app.GridLayout7 = uigridlayout(app.GridLayout4);
            app.GridLayout7.ColumnWidth = {'1x', '1x', '1x'};
            app.GridLayout7.RowHeight = {'1x'};
            app.GridLayout7.ColumnSpacing = 2;
            app.GridLayout7.RowSpacing = 1;
            app.GridLayout7.Padding = [0 0 0 0];
            app.GridLayout7.Layout.Row = 6;
            app.GridLayout7.Layout.Column = 1;

            % Create ButtonRun
            app.ButtonRun = uibutton(app.GridLayout7, 'push');
            app.ButtonRun.ButtonPushedFcn = createCallbackFcn(app, @ButtonRunPushed, true);
            app.ButtonRun.Layout.Row = 1;
            app.ButtonRun.Layout.Column = 2;
            app.ButtonRun.Text = 'Run';

            % Create ButtonGroup
            app.ButtonGroup = uibuttongroup(app.GridLayout7);
            app.ButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroupSelectionChanged, true);
            app.ButtonGroup.TitlePosition = 'centertop';
            app.ButtonGroup.Layout.Row = 1;
            app.ButtonGroup.Layout.Column = 1;

            % Create Button3_1
            app.Button3_1 = uiradiobutton(app.ButtonGroup);
            app.Button3_1.Text = '3:1';
            app.Button3_1.FontName = 'Arial';
            app.Button3_1.FontWeight = 'bold';
            app.Button3_1.Position = [0 20 58 22];
            app.Button3_1.Value = true;

            % Create Button1_3
            app.Button1_3 = uiradiobutton(app.ButtonGroup);
            app.Button1_3.Text = '1:3';
            app.Button1_3.FontName = 'Arial';
            app.Button1_3.FontWeight = 'bold';
            app.Button1_3.Position = [0 2 65 22];

            % Create ClearButton
            app.ClearButton = uibutton(app.GridLayout7, 'push');
            app.ClearButton.ButtonPushedFcn = createCallbackFcn(app, @ClearButtonPushed, true);
            app.ClearButton.Layout.Row = 1;
            app.ClearButton.Layout.Column = 3;
            app.ClearButton.Text = 'Clear';

            % Create Panel_3
            app.Panel_3 = uipanel(app.GridLayout4);
            app.Panel_3.TitlePosition = 'centertop';
            app.Panel_3.Layout.Row = 2;
            app.Panel_3.Layout.Column = 1;

            % Create MaximumPowerWEditFieldLabel
            app.MaximumPowerWEditFieldLabel = uilabel(app.Panel_3);
            app.MaximumPowerWEditFieldLabel.FontName = 'Arial';
            app.MaximumPowerWEditFieldLabel.FontWeight = 'bold';
            app.MaximumPowerWEditFieldLabel.Position = [8 11 119 22];
            app.MaximumPowerWEditFieldLabel.Text = 'Maximum Power(W)';

            % Create MaximumPowerWEditField
            app.MaximumPowerWEditField = uieditfield(app.Panel_3, 'numeric');
            app.MaximumPowerWEditField.Limits = [0 Inf];
            app.MaximumPowerWEditField.ValueChangedFcn = createCallbackFcn(app, @MaximumPowerWEditFieldValueChanged, true);
            app.MaximumPowerWEditField.Position = [126 11 46 22];

            % Create GridLayout5
            app.GridLayout5 = uigridlayout(app.GridLayout2);
            app.GridLayout5.ColumnWidth = {'1x'};
            app.GridLayout5.RowHeight = {'0.75x', '0.25x'};
            app.GridLayout5.ColumnSpacing = 0;
            app.GridLayout5.RowSpacing = 0;
            app.GridLayout5.Padding = [0 0 0 0];
            app.GridLayout5.Layout.Row = 1;
            app.GridLayout5.Layout.Column = 2;

            % Create GridLayout6
            app.GridLayout6 = uigridlayout(app.GridLayout5);
            app.GridLayout6.ColumnWidth = {'1x', '1x', '1x'};
            app.GridLayout6.RowHeight = {'1x'};
            app.GridLayout6.ColumnSpacing = 3;
            app.GridLayout6.RowSpacing = 3;
            app.GridLayout6.Padding = [3 3 3 3];
            app.GridLayout6.Layout.Row = 2;
            app.GridLayout6.Layout.Column = 1;

            % Create TreeSP
            app.TreeSP = uitree(app.GridLayout6, 'checkbox');
            app.TreeSP.Layout.Row = 1;
            app.TreeSP.Layout.Column = 1;

            % Create SeriesParallelNode
            app.SeriesParallelNode = uitreenode(app.TreeSP);
            app.SeriesParallelNode.Text = 'Series-Parallel';

            % Create VoltagesNode
            app.VoltagesNode = uitreenode(app.SeriesParallelNode);
            app.VoltagesNode.Text = 'Voltages';

            % Create SPVin
            app.SPVin = uitreenode(app.VoltagesNode);
            app.SPVin.Text = 'Vin';

            % Create SPVout
            app.SPVout = uitreenode(app.VoltagesNode);
            app.SPVout.Text = 'Vout';

            % Create CurrentsNode
            app.CurrentsNode = uitreenode(app.SeriesParallelNode);
            app.CurrentsNode.Text = 'Currents';

            % Create IinNode
            app.IinNode = uitreenode(app.CurrentsNode);
            app.IinNode.Text = 'Iin';

            % Create IoutNode
            app.IoutNode = uitreenode(app.CurrentsNode);
            app.IoutNode.Text = 'Iout';

            % Assign Checked Nodes
            app.TreeSP.CheckedNodesChangedFcn = createCallbackFcn(app, @TreeSPCheckedNodesChanged, true);

            % Create TreeLadder
            app.TreeLadder = uitree(app.GridLayout6, 'checkbox');
            app.TreeLadder.Layout.Row = 1;
            app.TreeLadder.Layout.Column = 2;

            % Create LadderNode
            app.LadderNode = uitreenode(app.TreeLadder);
            app.LadderNode.Text = 'Ladder';

            % Create VoltagesNode_2
            app.VoltagesNode_2 = uitreenode(app.LadderNode);
            app.VoltagesNode_2.Text = 'Voltages';

            % Create LadVin
            app.LadVin = uitreenode(app.VoltagesNode_2);
            app.LadVin.Text = 'Vin';

            % Create LadVout
            app.LadVout = uitreenode(app.VoltagesNode_2);
            app.LadVout.Text = 'Vout';

            % Create CurrentsNode_2
            app.CurrentsNode_2 = uitreenode(app.LadderNode);
            app.CurrentsNode_2.Text = 'Currents';

            % Create IinNode_2
            app.IinNode_2 = uitreenode(app.CurrentsNode_2);
            app.IinNode_2.Text = 'Iin';

            % Create IoutNode_2
            app.IoutNode_2 = uitreenode(app.CurrentsNode_2);
            app.IoutNode_2.Text = 'Iout';

            % Assign Checked Nodes
            app.TreeLadder.CheckedNodesChangedFcn = createCallbackFcn(app, @TreeLadderCheckedNodesChanged, true);

            % Create TreeDickson
            app.TreeDickson = uitree(app.GridLayout6, 'checkbox');
            app.TreeDickson.Layout.Row = 1;
            app.TreeDickson.Layout.Column = 3;

            % Create DicksonNode
            app.DicksonNode = uitreenode(app.TreeDickson);
            app.DicksonNode.Text = 'Dickson';

            % Create VoltajesNode
            app.VoltajesNode = uitreenode(app.DicksonNode);
            app.VoltajesNode.Text = 'Voltajes';

            % Create DckVin
            app.DckVin = uitreenode(app.VoltajesNode);
            app.DckVin.Text = 'Vin';

            % Create DckVout
            app.DckVout = uitreenode(app.VoltajesNode);
            app.DckVout.Text = 'Vout';

            % Create CurrentsNode_3
            app.CurrentsNode_3 = uitreenode(app.DicksonNode);
            app.CurrentsNode_3.Text = 'Currents';

            % Create IinNode_3
            app.IinNode_3 = uitreenode(app.CurrentsNode_3);
            app.IinNode_3.Text = 'Iin';

            % Create IoutNode_3
            app.IoutNode_3 = uitreenode(app.CurrentsNode_3);
            app.IoutNode_3.Text = 'Iout';

            % Assign Checked Nodes
            app.TreeDickson.CheckedNodesChangedFcn = createCallbackFcn(app, @TreeDicksonCheckedNodesChanged, true);

            % Create TabGroup
            app.TabGroup = uitabgroup(app.GridLayout5);
            app.TabGroup.Layout.Row = 1;
            app.TabGroup.Layout.Column = 1;

            % Create TimeResponseTab
            app.TimeResponseTab = uitab(app.TabGroup);
            app.TimeResponseTab.Title = 'Time Response';

            % Create UIAxes
            app.UIAxes = uiaxes(app.TimeResponseTab);
            title(app.UIAxes, 'Time Response - SC Converter')
            xlabel(app.UIAxes, 'Time(s)')
            ylabel(app.UIAxes, 'Voltage (V)')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.FontName = 'Arial';
            app.UIAxes.XLim = [0 Inf];
            app.UIAxes.GridLineStyle = '--';
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.Position = [0 0 442 207];

            % Create EfficiencyGraphTab
            app.EfficiencyGraphTab = uitab(app.TabGroup);
            app.EfficiencyGraphTab.Title = 'Efficiency Graph';

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.EfficiencyGraphTab);
            title(app.UIAxes2, 'Efficiency SCC Vs Frequency')
            xlabel(app.UIAxes2, 'Frequency(Hz)')
            ylabel(app.UIAxes2, 'Efficiency')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.FontName = 'Arial';
            app.UIAxes2.GridLineStyle = '--';
            app.UIAxes2.XGrid = 'on';
            app.UIAxes2.YGrid = 'on';
            app.UIAxes2.Position = [1 6 441 201];

            % Create PerformanceParametersTab
            app.PerformanceParametersTab = uitab(app.TabGroup);
            app.PerformanceParametersTab.Title = 'Performance Parameters';

            % Create GridLayout11
            app.GridLayout11 = uigridlayout(app.PerformanceParametersTab);
            app.GridLayout11.ColumnWidth = {'1x'};
            app.GridLayout11.ColumnSpacing = 0;
            app.GridLayout11.RowSpacing = 0;
            app.GridLayout11.Padding = [0 0 0 0];

            % Create Panel_7
            app.Panel_7 = uipanel(app.GridLayout11);
            app.Panel_7.TitlePosition = 'centertop';
            app.Panel_7.Title = '3:1';
            app.Panel_7.Layout.Row = 1;
            app.Panel_7.Layout.Column = 1;

            % Create GridLayout12
            app.GridLayout12 = uigridlayout(app.Panel_7);
            app.GridLayout12.ColumnWidth = {'1x'};
            app.GridLayout12.RowHeight = {'1x'};
            app.GridLayout12.ColumnSpacing = 0;
            app.GridLayout12.RowSpacing = 0;
            app.GridLayout12.Padding = [0 0 0 0];

            % Create UITable3_1
            app.UITable3_1 = uitable(app.GridLayout12);
            app.UITable3_1.ColumnName = {'Vout(V)'; 'R_line(ΔVout/ΔVin)'; 'R_load(ΔVout/ΔIoad)'; 'P_loos(W)'; 'Efficiency'; 'Mssl'; 'Mfsl'};
            app.UITable3_1.RowName = {};
            app.UITable3_1.FontName = 'Arial';
            app.UITable3_1.Layout.Row = 1;
            app.UITable3_1.Layout.Column = 1;

            % Create Panel_8
            app.Panel_8 = uipanel(app.GridLayout11);
            app.Panel_8.TitlePosition = 'centertop';
            app.Panel_8.Title = '1:3';
            app.Panel_8.Layout.Row = 2;
            app.Panel_8.Layout.Column = 1;

            % Create GridLayout13
            app.GridLayout13 = uigridlayout(app.Panel_8);
            app.GridLayout13.ColumnWidth = {'1x'};
            app.GridLayout13.RowHeight = {'1x'};
            app.GridLayout13.ColumnSpacing = 0;
            app.GridLayout13.RowSpacing = 0;
            app.GridLayout13.Padding = [0 0 0 0];

            % Create UITable1_3
            app.UITable1_3 = uitable(app.GridLayout13);
            app.UITable1_3.ColumnName = {'Vout(V)'; 'R_line(ΔVout/ΔVin)'; 'R_load(ΔVout/ΔIoad)'; 'P_loos(W)'; 'Efficiency'; 'Mssl'; 'Mfsl'};
            app.UITable1_3.RowName = {};
            app.UITable1_3.Layout.Row = 1;
            app.UITable1_3.Layout.Column = 1;

            % Create PanelTop
            app.PanelTop = uipanel(app.GridLayout);
            app.PanelTop.TitlePosition = 'centertop';
            app.PanelTop.Layout.Row = 1;
            app.PanelTop.Layout.Column = 1;

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.PanelTop);
            app.GridLayout3.ColumnWidth = {'0.1x', '1x', '0.1x'};
            app.GridLayout3.RowHeight = {'1x'};

            % Create BidirectionalSwitchedCapacitorConverterAnalyzerLabel
            app.BidirectionalSwitchedCapacitorConverterAnalyzerLabel = uilabel(app.GridLayout3);
            app.BidirectionalSwitchedCapacitorConverterAnalyzerLabel.HorizontalAlignment = 'center';
            app.BidirectionalSwitchedCapacitorConverterAnalyzerLabel.FontName = 'Arial Black';
            app.BidirectionalSwitchedCapacitorConverterAnalyzerLabel.FontSize = 18;
            app.BidirectionalSwitchedCapacitorConverterAnalyzerLabel.FontWeight = 'bold';
            app.BidirectionalSwitchedCapacitorConverterAnalyzerLabel.Layout.Row = 1;
            app.BidirectionalSwitchedCapacitorConverterAnalyzerLabel.Layout.Column = 2;
            app.BidirectionalSwitchedCapacitorConverterAnalyzerLabel.Interpreter = 'tex';
            app.BidirectionalSwitchedCapacitorConverterAnalyzerLabel.Text = {'Bidirectional Switched-Capacitor Converter Analyzer'; '3:1 ↔ 1:3 Topology Comparison Tool'};

            % Create GridLayout9
            app.GridLayout9 = uigridlayout(app.GridLayout3);
            app.GridLayout9.ColumnSpacing = 0;
            app.GridLayout9.RowSpacing = 0;
            app.GridLayout9.Padding = [1 1 1 1];
            app.GridLayout9.Layout.Row = 1;
            app.GridLayout9.Layout.Column = 3;

            % Create Lamp
            app.Lamp = uilamp(app.GridLayout9);
            app.Lamp.Layout.Row = 2;
            app.Lamp.Layout.Column = 2;

            % Create GridLayout10
            app.GridLayout10 = uigridlayout(app.GridLayout);
            app.GridLayout10.ColumnWidth = {'1x'};
            app.GridLayout10.RowHeight = {'1x'};
            app.GridLayout10.ColumnSpacing = 0;
            app.GridLayout10.RowSpacing = 0;
            app.GridLayout10.Padding = [0 0 0 0];
            app.GridLayout10.Layout.Row = 3;
            app.GridLayout10.Layout.Column = 1;

            % Create PanelBottom
            app.PanelBottom = uipanel(app.GridLayout10);
            app.PanelBottom.TitlePosition = 'centertop';
            app.PanelBottom.Title = 'Optimized Components';
            app.PanelBottom.Layout.Row = 1;
            app.PanelBottom.Layout.Column = 1;
            app.PanelBottom.FontName = 'Arial';
            app.PanelBottom.FontWeight = 'bold';

            % Create GridLayout8
            app.GridLayout8 = uigridlayout(app.PanelBottom);
            app.GridLayout8.ColumnWidth = {'0.4x', '0.6x'};
            app.GridLayout8.RowHeight = {'1x'};
            app.GridLayout8.ColumnSpacing = 1;
            app.GridLayout8.RowSpacing = 1;
            app.GridLayout8.Padding = [1 1 1 1];

            % Create UITableCapacitor
            app.UITableCapacitor = uitable(app.GridLayout8);
            app.UITableCapacitor.ColumnName = {'C1'; 'C2'; 'C3'; '3:1Cout'; '1:3Cout'};
            app.UITableCapacitor.RowName = {};
            app.UITableCapacitor.SelectionType = 'row';
            app.UITableCapacitor.Layout.Row = 1;
            app.UITableCapacitor.Layout.Column = 1;

            % Create UITableSwitch
            app.UITableSwitch = uitable(app.GridLayout8);
            app.UITableSwitch.ColumnName = {'S1'; 'S2'; 'S3'; 'S4'; 'S5'; 'S6'; 'S7'};
            app.UITableSwitch.RowName = {};
            app.UITableSwitch.Layout.Row = 1;
            app.UITableSwitch.Layout.Column = 2;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SC_Comparador_final

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
