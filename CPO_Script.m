% ***********************************************************************************
%         C A R D I A C   P O W E R   O U T P U T   C A L C U L A T I O N
% ***********************************************************************************
%
%   This script uses a patient specific cardiovascular systems model which is
%   a simplification of the Smith el al. model (Med Eng Phys 26:131, 2004) to 
%   calculate the cardiac power output from the model optimized to the set of 
%   longitudinal data for a given patient.
%
%   Model originally created on     29 March 2019
%   Model last modfied on           3  April 2019

%   Developed by        Brian Carlson
%                       Physiological Systems Dynamics Laboratory
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
%
% ***********************************************************************************
%  Start of                  C P O   C A L C U L A T I O N    S C R I P T
% ***********************************************************************************

%% **********************************************************************************
%  Run Calc for             C P O   C A L C U L A T I O N    S C R I P T
% ***********************************************************************************
    

    Patient_Num = [266 363 456 558 572];
    PlotSimFits_Flag = 0;
    VisCPOCalcs_Flag = 0;
    
    Conv1 = 2.22204e-6;                                 % From (mmHg*mL)/min to W
    
    Num_Pats = size(Patient_Num,2);
    Num_RHCs = zeros(1,Num_Pats);
    for i = 1:Num_Pats
        
        % Get RHC data and optimized parameters for patient i
        Lngtdnl_OptimP{i} = ...
            importdata(['LngtdnlOptP_Patient', ...
            num2str(Patient_Num(i)), '.txt']);
        
        % Number of RHCs in the longitudinal dataset for patient i
        Num_RHCs(i) = size(Lngtdnl_OptimP{i}.data,2);
        
        % Step through and simulate each RHC dataset in the logitudinal dataset
        for j = 1:Num_RHCs(i)
        
            % RHC PARAMETERS
            Days_PTx = Lngtdnl_OptimP{i}.data(1,j);    % Days post transplant 
            HR = Lngtdnl_OptimP{i}.data(2,j);          % Heart rate (beats/min)
            Hgt = Lngtdnl_OptimP{i}.data(3,j);         % Patient height (cm)
            Wght = Lngtdnl_OptimP{i}.data(4,j);        % Patient weight (kg)
            Gender = Lngtdnl_OptimP{i}.data(5,j);      % Ptnt gndr (1-male, 2-female)
            % MODEL PARAMETERS
            % Left ventricle free wall parameters
            E_es_lvf = Lngtdnl_OptimP{i}.data(6,j);    % LV free wl elstnce (kPa/mL)
            V_d_lvf = Lngtdnl_OptimP{i}.data(7,j);     % LV ES zero P volume (mL)
            P_0_lvf = Lngtdnl_OptimP{i}.data(8,j);     % LV ED pressure param (kPa)
            lambda_lvf = Lngtdnl_OptimP{i}.data(9,j);  % LV ED pressure param (1/mL)
            V_0_lvf = Lngtdnl_OptimP{i}.data(10,j);    % LV ED pressure param (mL)
            % Right ventricle free wall parameters
            E_es_rvf = Lngtdnl_OptimP{i}.data(11,j);   % RV free wl elstnce (kPa/mL)
            V_d_rvf = Lngtdnl_OptimP{i}.data(12,j);    % RV ES zero P volume (mL)
            P_0_rvf = Lngtdnl_OptimP{i}.data(13,j);    % RV ED pressure param (kPa)
            lambda_rvf = Lngtdnl_OptimP{i}.data(14,j); % RV ED pressure param (1/mL)
            V_0_rvf = Lngtdnl_OptimP{i}.data(15,j);    % RV ED pressure param (mL)
            % Pulmonary artery and vein parameters
            E_es_pa = Lngtdnl_OptimP{i}.data(16,j);    % Pulm art elstnce (kPa/mL)
            E_es_pu = Lngtdnl_OptimP{i}.data(17,j);    % Pulm vein elstnce (kPa/mL)
            R_pul = Lngtdnl_OptimP{i}.data(18,j);      % Pulm vasc rsstnce (kPa*s/mL)
            P_th = Lngtdnl_OptimP{i}.data(19,j);       % Mean thrcc pressure (mmHg)
            % Aortic and vena cava parameters
            E_es_ao = Lngtdnl_OptimP{i}.data(20,j);    % Aorta elastance (kPa/mL)
            E_es_vc = Lngtdnl_OptimP{i}.data(21,j);    % Vena cava elastance (kPa/mL)
            R_sys = Lngtdnl_OptimP{i}.data(22,j);      % Syst art rsstnce (kPa*s/mL)
            % Heart valve parameters
            R_mt = Lngtdnl_OptimP{i}.data(23,j);       % Mtrl vlv resist (kPa*s/mL)
            R_av = Lngtdnl_OptimP{i}.data(24,j);       % Artc valve resist (kPa*s/mL)
            R_tc = Lngtdnl_OptimP{i}.data(25,j);       % Trcspd vlv resist (kPa*s/mL)
            R_pv = Lngtdnl_OptimP{i}.data(26,j);       % Pulm vlv resist (kPa*s/mL)
            
            RHCParam_Values = {Days_PTx HR Hgt Wght Gender};
            RHCParam_Fields = {'Days_PTx' 'HR' 'Hgt' 'Wght' 'Gender'};
            RHCParam_Struct = cell2struct(RHCParam_Values, ...
                RHCParam_Fields,2);
            
            CVParam_Values = {E_es_lvf V_d_lvf P_0_lvf lambda_lvf V_0_lvf ...
                E_es_rvf V_d_rvf P_0_rvf lambda_rvf V_0_rvf E_es_pa E_es_pu ...
                R_pul P_th E_es_ao E_es_vc R_sys R_mt R_av R_tc R_pv};
            CVParam_Fields = {'E_es_lvf' 'V_d_lvf' 'P_0_lvf' 'lambda_lvf' ...
                'V_0_lvf' 'E_es_rvf' 'V_d_rvf' 'P_0_rvf' 'lambda_rvf' ...
                'V_0_rvf' 'E_es_pa' 'E_es_pu' 'R_pul' 'P_th' 'E_es_ao' ...
                'E_es_vc' 'R_sys' 'R_mt' 'R_av' 'R_tc' 'R_pv'};
            CVParam_Struct = cell2struct(CVParam_Values, ...
                CVParam_Fields,2);
            
            % TOTAL AND CIRCULATING BLOOD VOLUME
            % Calculate total blood volume based on height (in cm), 
            %  weight (in kg) and gender (1-male, 2-female). This 
            %  expression is from Nadler et al. Surgery 51:224,1962.
            if (Gender == 1)
                TotBV = ((0.3669 * (Hgt/100)^3) + (0.03219 * Wght) + 0.6041) * 1000;
            else
                TotBV = ((0.3561 * (Hgt/100)^3) + (0.03308 * Wght) + 0.1833) * 1000;
            end
            
            % The original Smith model only circulated a portion of the blood 
            %  so aortic pressure dynamics are not lumped into a general arterial 
            %  systemic compartment. Assuming they were simulating a typical 5000 mL
            %  total blood volume they included only 1500 mL (or 30%) in the 
            %  circulating volume therefore we will multiply our calculated TotBV 
            %  value by 0.3 to yield circulating blood volume. 
            CircBV = 0.30 * TotBV; 
            
            % INITIAL CONDITIONS
            %  Note that initial volume division is scaled as a fraction
            %  of circulating blood volume calculated earlier
            V_lv0 = (94.6812/1500) * CircBV;
            V_rv0 = (90.7302/1500) * CircBV;
            V_pa0 = (43.0123/1500) * CircBV;
            V_pu0 = (808.458/1500) * CircBV;
            V_ao0 = (133.338/1500) * CircBV;
            V_vc0 = (329.780/1500) * CircBV;
            
            X0(1) = V_lv0;
            X0(2) = V_rv0;
            X0(3) = V_pa0;
            X0(4) = V_pu0;
            X0(5) = V_ao0;
            X0(6) = V_vc0;
            
            % SIMULATION CONDITIONS
            % Simulation time scales
            NumBeats_SS = 20;
            NumBeats_Sim = 5; 
            Per = 60 / HR;                            % Period (sec/beat)
            TSpan_SS = [0 NumBeats_SS * Per];
            TSpan_Sim = [0 NumBeats_Sim * Per];
            
            % MODEL SIMULATION
            % Set ODE/DAE options
            ODE_Opts = odeset('AbsTol',1e-7,'RelTol',1e-4); 
            % Solve over the steady state time span with ode15s
            [T_Sim_SS,X_Sim_SS] = ...
                ode15s(@dXdT_HTxCV,TSpan_SS,X0, ...
                ODE_Opts,RHCParam_Struct,CVParam_Struct);
            % Now solve over the simulation time span
            [T_Sim{i,j},X_Sim{i,j}] = ...
                ode15s(@dXdT_HTxCV,TSpan_Sim, ...
                X_Sim_SS(end,:),ODE_Opts, ...
                RHCParam_Struct,CVParam_Struct);
            
            % EXTRACTING PRESSURES FROM MODEL SOLUTION 
            Num_TSim = size(T_Sim{i,j},1);          % Number of timepoints
            % Calculate and store intermediate LV and RV pressures
            for k = 1:Num_TSim
                VarOut = dXdT_HTxCV(T_Sim{i,j}(k), ...
                    X_Sim{i,j}(k,:),RHCParam_Struct, ...
                    CVParam_Struct,1);
                P_LVSim{i,j}(k) = VarOut(1);
                P_RVSim{i,j}(k) = VarOut(2);
                P_AOSim{i,j}(k) = VarOut(3);
                P_VCSim{i,j}(k) = VarOut(4);
                P_PASim{i,j}(k) = VarOut(5);
                P_PUSim{i,j}(k) = VarOut(6);
            end
            
        end
        
    end
    
    
%% **********************************************************************************
%  Power Output Calc for    C P O   C A L C U L A T I O N    S C R I P T
% ***********************************************************************************
    
    for i = 1:Num_Pats
        
        for j = 1:Num_RHCs(i)
            
            % Find max and min left and right ventricular 
            %  pressures for each beat
            P_LVMaxMinInd = [];
            P_RVMaxMinInd = [];
            Num_PLVPts = size(P_LVSim{i,j},2);
            P_LVSimMax = max(P_LVSim{i,j});
            P_LVSimMin = min(P_LVSim{i,j});
            P_LVSimMean = (P_LVSimMax + P_LVSimMin) / 2;
            P_RVSimMax = max(P_RVSim{i,j});
            P_RVSimMin = min(P_RVSim{i,j});
            P_RVSimMean = (P_RVSimMax + P_RVSimMin) / 2;
            P_LVSimDiff = diff(P_LVSim{i,j});                   % Diffs of P to find
            P_RVSimDiff = diff(P_RVSim{i,j});                   %  change in sign
            for k = 2:Num_PLVPts-1
                if (size(P_LVMaxMinInd) == 0)
                    SkipLVMaxMin_Flag = 0;  
                else
                    P_LVMaxDelta = abs(P_LVSimMax -  ...
                        P_LVSim{i,j}(P_LVMaxMinInd(end)));
                    P_LVMinDelta = abs(P_LVSim{i,j} ...
                        (P_LVMaxMinInd(end)) - P_LVSimMin);
                    if (P_LVMaxDelta >= P_LVMinDelta)           % Search for max
                        if (P_LVSim{i,j}(k) >= P_LVSimMean)
                            SkipLVMaxMin_Flag = 0;
                        else
                            SkipLVMaxMin_Flag = 1;
                        end
                    else                                        % Search for min
                        if (P_LVSim{i,j}(k) <= P_LVSimMean)
                            SkipLVMaxMin_Flag = 0;
                        else
                            SkipLVMaxMin_Flag = 1;
                        end
                    end
                end
                
                if (size(P_RVMaxMinInd) == 0)
                    SkipRVMaxMin_Flag = 0;  
                else    
                    P_RVMaxDelta = abs(P_RVSimMax - ...
                        P_RVSim{i,j}(P_RVMaxMinInd(end)));
                    P_RVMinDelta = abs(P_RVSim{i,j} ...
                        (P_RVMaxMinInd(end)) - P_RVSimMin);
                    if (P_RVMaxDelta >= P_RVMinDelta)           % Search for max
                        if (P_RVSim{i,j}(k) >= P_RVSimMean)
                            SkipRVMaxMin_Flag = 0;
                        else
                            SkipRVMaxMin_Flag = 1;
                        end
                    else                                        % Search for min
                        if (P_RVSim{i,j}(k) <= P_RVSimMean)
                            SkipRVMaxMin_Flag = 0;
                        else
                            SkipRVMaxMin_Flag = 1;
                        end    
                    end
                end
                    
                if ((P_LVSimDiff(k-1)*P_LVSimDiff(k)) <= 0 ...  % Store indices
                        && SkipLVMaxMin_Flag == 0)              %  where sign
                    P_LVMaxMinInd = [P_LVMaxMinInd k];          %  changes for
                end                                             %  both LV and RV
                if ((P_RVSimDiff(k-1)*P_RVSimDiff(k)) <= 0 ...  %  as long as we 
                        && SkipRVMaxMin_Flag == 0)              %  are far enough 
                    P_RVMaxMinInd = [P_RVMaxMinInd k];          %  away from last 
                end                                             %  min/max
            end
            
            % Now use the max and min pressures to form brackets to search 
            %  for the max and min volumes for both LV and RV
                     
            % Find Left Ventricle max and min volumes
            %  and the indices associated with them
            P_LVMean = mean(P_LVSim{i,j});                      % Mean P_LV           
            Num_PLVMaxMin = size(P_LVMaxMinInd,2);              % Num LVP max,mins    
            if (P_LVSim{i,j}(P_LVMaxMinInd(1)) >= P_LVMean)     % Check to see
                Min_Flag = 1;                                   %  if the first LVP
            else                                                %  max/min index is
                Min_Flag = 0;                                   %  a max or a min
            end
            V_LVLowInd = [];                                    % Init max/min V_LV 
            V_LVHiInd = [];                                     %  index vectors
            for k = 1:Num_PLVMaxMin-1                           % Search between
                V_LVBrackLowInd = P_LVMaxMinInd(k);             %  successive
                V_LVBrackHiInd = P_LVMaxMinInd(k+1);            %  max/min P_LV
                if (Min_Flag == 1)                              %  indicies to find
                    V_LVLowVal = min(X_Sim{i,j} ...             %  max/min V_LVs
                        (V_LVBrackLowInd:V_LVBrackHiInd,1));    
                    V_LVLowInd = ...                            % Here we find the
                        [V_LVLowInd (V_LVBrackLowInd + ...      %  first minimum  
                        find(X_Sim{i,j}(V_LVBrackLowInd: ...    %  V_LV and store
                        V_LVBrackHiInd,1) == ...                %  it's index
                        V_LVLowVal,1,'first'))];
                    Min_Flag = 0;                               % Look for max next
                else
                    V_LVHiVal = max(X_Sim{i,j} ...  
                        (V_LVBrackLowInd:V_LVBrackHiInd,1));
                    V_LVHiInd =  ...                            % Here we find the
                        [V_LVHiInd (V_LVBrackLowInd + ...       %  first maximum
                        find(X_Sim{i,j}(V_LVBrackLowInd: ...    %  V_LV and store
                        V_LVBrackHiInd,1) == ...                %  it's index
                        V_LVHiVal,1,'first'))];
                    Min_Flag = 1;                               % Look for min next
                end
            end
            % Find Right Ventricle max and min volumes
            %  and the indices associated with them   
            P_RVMean = mean(P_RVSim{i,j});                      % Mean P_RV
            Num_PRVMaxMin = size(P_RVMaxMinInd,2);              % Num RVP max,mins
            if (P_RVSim{i,j}(P_RVMaxMinInd(1)) >= P_RVMean)     % Check to see
                Min_Flag = 1;                                   %  if the first RVP
            else                                                %  max/min index is
                Min_Flag = 0;                                   %  a max or a min
            end
            V_RVLowInd = [];                                    % Init max/min V_RV
            V_RVHiInd = [];                                     %  index vectors
            for k = 1:Num_PRVMaxMin-1                           % Search between 
                V_RVBrackLowInd = P_RVMaxMinInd(k);             %  successive 
                V_RVBrackHiInd = P_RVMaxMinInd(k+1);            %  max/min P_RV
                if (Min_Flag == 1)                              %  indices to find
                    V_RVLowVal = min(X_Sim{i,j} ...             %  max/min V_RVs
                        (V_RVBrackLowInd:V_RVBrackHiInd,2));
                    V_RVLowInd = ...                            % Here we find the
                        [V_RVLowInd (V_RVBrackLowInd + ...      %  first minimum
                        find(X_Sim{i,j}(V_RVBrackLowInd: ...    %  V_RV and store
                        V_RVBrackHiInd,2) == ...                %  it's index
                        V_RVLowVal,1,'first'))];
                    Min_Flag = 0;                               % Look for max next
                else
                    V_RVHiVal = max(X_Sim{i,j} ...
                        (V_RVBrackLowInd:V_RVBrackHiInd,2));
                    V_RVHiInd =  ...                            % Here we find the
                        [V_RVHiInd (V_RVBrackLowInd + ...       %  first maximum
                        find(X_Sim{i,j}(V_RVBrackLowInd: ...    %  V_LV and store
                        V_RVBrackHiInd,2) == ...                %  it's index
                        V_RVHiVal,1,'first'))];
                    Min_Flag = 1;                               % Look for max next
                end
            end
            
            % Now find calculate the power output over each 
            %  cycle for both the left and right ventricle
            HR = Lngtdnl_OptimP{i}.data(2,j);                   % Hrt rt (beats/min)
            Num_CPOBeats = (Num_PLVMaxMin/2) - 2;               % Num of beats
            if (V_LVLowInd(1) <= V_LVHiInd(1))                  % Make sure to start
                V_LowIndShift = 1;                              %  on contraction
            else                                                %  which is the same
                V_LowIndShift = 0;                              %  for both LV and
            end                                                 %  RV
            
            for k = 1:Num_CPOBeats
                % Left ventricular CPO calculation
                LV_PVAContStrtInd = V_LVHiInd(k);               % Set cont start 
                LV_PVAContEndInd = ...                          %  and end indices
                    V_LVLowInd(k+V_LowIndShift);
                
                LV_PVACont = abs(trapz(X_Sim{i,j} ...           % Integrate the 
                    (LV_PVAContStrtInd: ...                     %  contraction
                    LV_PVAContEndInd,1), ...                    %  pressure volume
                    P_LVSim{i,j}(LV_PVAContStrtInd: ...         %  area
                    LV_PVAContEndInd)));
                
                LV_PVAFillStrtInd = ...                         % Set fill start
                    V_LVLowInd(k+V_LowIndShift);                %  and end indices
                LV_PVAFillEndInd = V_LVHiInd(k+1);     
                LV_PVAFill = abs(trapz(X_Sim{i,j} ...           % Integrate the 
                    (LV_PVAFillStrtInd: ...                     %  filling pressure
                    LV_PVAFillEndInd,1), ...                    %  volume area
                    P_LVSim{i,j}(LV_PVAFillStrtInd: ...
                    LV_PVAFillEndInd)));
                
                LV_PVA = LV_PVACont - LV_PVAFill;               % Subtract fill and
                LV_CPO{i,j}(k) = LV_PVA * HR * Conv1;           %  mult by HR (Watts)
                
                % Right ventricular CPO calculation
                RV_PVAContStrtInd = V_RVHiInd(k);               % Set cont start 
                RV_PVAContEndInd = ...                          %  and end indices
                    V_RVLowInd(k+V_LowIndShift);
                RV_PVACont = abs(trapz(X_Sim{i,j} ...           % Integrate the 
                    (RV_PVAContStrtInd: ...                     %  contraction
                    RV_PVAContEndInd,2), ...                    %  pressure volume
                    P_RVSim{i,j}(RV_PVAContStrtInd: ...         %  area
                    RV_PVAContEndInd)));
                
                RV_PVAFillStrtInd = ...                         % Set fill start
                    V_RVLowInd(k+V_LowIndShift);                %  and end indices
                RV_PVAFillEndInd = V_RVHiInd(k+1);     
                RV_PVAFill = abs(trapz(X_Sim{i,j} ...           % Integrate the 
                    (RV_PVAFillStrtInd: ...                     %  filling pressure
                    RV_PVAFillEndInd,2), ...                    %  volume area
                    P_RVSim{i,j}(RV_PVAFillStrtInd: ...
                    RV_PVAFillEndInd)));
                
                RV_PVA = RV_PVACont - RV_PVAFill;               % Subtract fill and
                RV_CPO{i,j}(k) = RV_PVA * HR * Conv1;           %  mult by HR (Watts)
                
            end
            
            if (VisCPOCalcs_Flag == 1)
                
                figure(3)
                subplot(2,1,1)
                plot(T_Sim{i,j},P_LVSim{i,j},'-b')
                hold on
                for k = 1:Num_PLVMaxMin-1
                    plot(T_Sim{i,j}(P_LVMaxMinInd(k)), ...
                        P_LVSim{i,j}(P_LVMaxMinInd(k)),'*r')
                end
                plot(T_Sim{i,j},X_Sim{i,j}(:,1),'-k')
                for k = 1:Num_PLVMaxMin/2
                    plot(T_Sim{i,j}(V_LVLowInd(k)), ...
                        X_Sim{i,j}(V_LVLowInd(k),1),'+r')
                end
                for k = 1:(Num_PLVMaxMin/2)-1
                    plot(T_Sim{i,j}(V_LVHiInd(k)), ...
                        X_Sim{i,j}(V_LVHiInd(k),1),'or')
                end
            
                subplot(2,1,2)
                plot(T_Sim{i,j},P_RVSim{i,j},'-g')
                hold on
                for k = 1:Num_PRVMaxMin
                    plot(T_Sim{i,j}(P_RVMaxMinInd(k)), ...
                        P_RVSim{i,j}(P_RVMaxMinInd(k)),'*r')
                end
                plot(T_Sim{i,j},X_Sim{i,j}(:,2),'-m')
                for k = 1:Num_PRVMaxMin/2
                    plot(T_Sim{i,j}(V_RVLowInd(k)), ...
                        X_Sim{i,j}(V_RVLowInd(k),2),'+r')
                end
                for k = 1:(Num_PRVMaxMin/2)-1
                    plot(T_Sim{i,j}(V_RVHiInd(k)), ...
                        X_Sim{i,j}(V_RVHiInd(k),2),'or')
                end
            
            end
            
        end
        
    end
            
    
    
%% **********************************************************************************
%  Plot Figures for         C P O   C A L C U L A T I O N    S C R I P T
% ***********************************************************************************
    
    Fig_Num = 0;
    ScrSize = get(0,'ScreenSize');              % Getting screen size
    
    % PLOTTING SIMULATION FIT TO EACH DATA SET
    if (PlotSimFits_Flag == 1)
        for i = 1:Num_Pats
        
            % Get RHC data and optimized parameters for patient i
            Lngtdnl_RHC = ...
                importdata(['LngtdnlRHCData_Patient', ...
                num2str(Patient_Num(i)), '.txt']);
        
            for j = 1:Num_RHCs(i)

                P_RVsyst = Lngtdnl_RHC.data(1,j);       % RV systlc pressure (mmHg)
                P_RVdiast = Lngtdnl_RHC.data(2,j);      % RV diastlc pressure (mmHg)
                P_PAsyst = Lngtdnl_RHC.data(3,j);       % Pulm art systlc prsr (mmHg)
                P_PAdiast = Lngtdnl_RHC.data(4,j);      % Pulm art dstlc prsr (mmHg)
                P_PCWave = Lngtdnl_RHC.data(5,j);       % Pulm cap wdg prsr (mmHg)
                P_AOsyst = Lngtdnl_RHC.data(6,j);       % Aortic systlc prsr (mmHg)
                P_AOdiast = Lngtdnl_RHC.data(7,j);      % Aortic dstlc prsr (mmHg)
                CO_Fick = Lngtdnl_RHC.data(8,j);        % Crdc otpt by Fick (L/min)

                Fig_Num = Fig_Num + 1;

                figure(Fig_Num)
                % Right ventricular pressure
                plot(T_Sim{i,j},P_RVSim{i,j},'-g', ...      % P_RV simulation
                    'LineWidth',3, ...
                    'DisplayName','P_{RV}')
                hold on
                plot([T_Sim{i,j}(1),T_Sim{i,j}(end)], ...   % P_RV systole data
                    [P_RVsyst,P_RVsyst],'-.g', ...
                    'LineWidth',1.5, ...
                    'HandleVisibility','off')
                plot([T_Sim{i,j}(1),T_Sim{i,j}(end)], ...   % P_RV diastole data
                    [P_RVdiast,P_RVdiast],':g', ...
                    'LineWidth',1.5, ...
                    'HandleVisibility','off')
                % Aortic pressure
                plot(T_Sim{i,j},P_AOSim{i,j},'-r', ...      % P_AO simulation
                    'LineWidth',3, ...
                    'DisplayName','P_{AO}')
                plot([T_Sim{i,j}(1),T_Sim{i,j}(end)], ...   % P_AO systole data
                    [P_AOsyst,P_AOsyst],'-.r', ...
                    'LineWidth',1.5, ...
                    'HandleVisibility','off')
                plot([T_Sim{i,j}(1),T_Sim{i,j}(end)], ...   % P_AO diastole data
                    [P_AOdiast,P_AOdiast],':r', ...
                    'LineWidth',1.5, ...
                    'HandleVisibility','off')
                % Pulmonary artery pressure
                plot(T_Sim{i,j},P_PASim{i,j},'-b', ...      % P_PA simulation
                    'LineWidth',3, ...
                    'DisplayName','P_{PA}')
                plot([T_Sim{i,j}(1),T_Sim{i,j}(end)], ...   % P_PA systole data
                    [P_PAsyst,P_PAsyst],'-.b', ...
                    'LineWidth',1.5, ...
                    'HandleVisibility','off')
                plot([T_Sim{i,j}(1),T_Sim{i,j}(end)], ...   % P_PA diastole data
                    [P_PAdiast,P_PAdiast],':b', ...
                    'LineWidth',1.5, ...
                    'HandleVisibility','off')
                % Pulmonary vein pressure
                plot(T_Sim{i,j},P_PUSim{i,j},'-c', ...      % P_PU simulation
                    'LineWidth',3, ...
                    'DisplayName','P_{PU}')
                plot([T_Sim{i,j}(1),T_Sim{i,j}(end)], ...   % P_PCW average data
                    [P_PCWave,P_PCWave],'-.c', ...
                    'LineWidth',1.5, ...
                    'HandleVisibility','off')
                % Cardiac output 
                PMax_Figi = 1.45 * ...
                    max(P_AOsyst,max(P_AOSim{i,j}));
    %             text(0.5,0.95*PMax_SP1, ...                 % CO data
    %                 ['CO Data = ' num2str(CO_Fick)])
    %             text(0.5,0.90*PMax_SP1, ...
    %                 ['CO Sim  = ' num2str(CO_Sim)])         % CO simulation

                xlim([0 5])                                     % Format figure
                ylim([-20 PMax_Figi])
                LegHndl_Figi = legend('show');
                set(LegHndl_Figi,'Box','off','FontSize',8)
                set(gca,'FontSize',14,'FontWeight', ...
                    'bold','Box','off')
                xlabel('Time (sec)','FontSize',20, ...
                    'FontWeight','bold')
                ylabel('Pressures (mmHg)', ...
                    'FontSize',20,'FontWeight','bold')

            end

        end
        
    end
            
    % PLOT LV AND RV PRESSURE VOLUME LOOPS FOR EACH PATIENT  
    LnColor_Str = {'m' 'r' 'y' 'g' 'c' 'b' 'k' 'm'};
    for i = 1:Num_Pats
        
        Fig_Num = Fig_Num + 1;
        PV_Fig(Fig_Num) = figure('Position', ...                % Positioning the
        [ScrSize(3)/50 ScrSize(4)/50 ...                        %  figure on the 
        ScrSize(3)/3.25 ScrSize(4)/1.10]);                         %  screen
    
        PVLSupTtl_Hndl = ...
            suptitle(['Patient Number ' num2str(Patient_Num(i))]);
        set(PVLSupTtl_Hndl,'FontSize',20,'FontWeight','bold')
        
        LVVol_Max = 0;
        RVVol_Max = 0;
        PLV_Max = 0;
        Lgnd_Strgs = cell(1,Num_RHCs(i));
        
        % LV PV Loops
        subplot(2,1,1)
        for j = 1:Num_RHCs(i)

            LVVol_Max = ...
                max(max(X_Sim{i,j}(:,1)),LVVol_Max);
            RVVol_Max = ...
                max(max(X_Sim{i,j}(:,2)),RVVol_Max);
            PLV_Max = ...
                max(max(P_LVSim{i,j}),PLV_Max);
            PLV_Min = -10;

            plot(X_Sim{i,j}(:,1), ...                           % V_LV vs P_LV simul
                P_LVSim{i,j}, ...
                ['-' LnColor_Str{j}], ...
                'LineWidth',3)
            hold on
          
            Days_PTx = Lngtdnl_OptimP{i}.data(1,j);
            Lgnd_Strgs{j} = num2str(Days_PTx);

        end
        
        %xlim([0 LVVol_Max*1.20])                                % Format subplot
        xlim([50 225])
        ylim([PLV_Min 1.15*PLV_Max])
        LegHndlFigi_SP1 = legend(Lgnd_Strgs{:});
        set(LegHndlFigi_SP1,'Box','off','FontSize',8)
        set(gca,'FontSize',14, ...
            'FontWeight','bold','Box','off')
        xlabel('LV Volume (mL)', ...
            'FontSize',18,'FontWeight','bold')
        ylabel('LV Pressure (mmHg)', ...
            'FontSize',18,'FontWeight','bold')
        
         % RV PV Loops
        subplot(2,1,2)
        for j = 1:Num_RHCs(i)
            LVVol_Max = max(max(X_Sim{i,j}(:,1)),LVVol_Max);
            RVVol_Max = max(max(X_Sim{i,j}(:,2)),RVVol_Max);
            PLV_Max = max(max(P_LVSim{i,j}),PLV_Max);
            PLV_Min = -10;
            
            plot(X_Sim{i,j}(:,2), ...                           % V_RV vs P_RV simul
                P_RVSim{i,j}, ...
                ['-' LnColor_Str{j}], ...
                'LineWidth',3)
            hold on
        end
        
        %xlim([0 LVVol_Max*1.20])                                % Format subplot
        xlim([50 225])
        ylim([PLV_Min 1.15*PLV_Max])
        LegHndlFigi_SP2 = legend(Lgnd_Strgs{:});
        set(LegHndlFigi_SP2,'Box','off','FontSize',8)
        set(gca,'FontSize',14, ...
            'FontWeight','bold','Box','off')
        xlabel('RV Volume (mL)', ...
            'FontSize',18,'FontWeight','bold')
        ylabel('RV Pressure (mmHg)', ...
            'FontSize',18,'FontWeight','bold')
        
    end
    
    % PLOT TREND IN LV AND RV POWER OUTPUT FOR EACH PATIENT
    for i = 1:Num_Pats
        
        Fig_Num = Fig_Num + 1;
        CPO_Fig(Fig_Num) = figure('Position', ...               % Positioning the
        [ScrSize(3)/50 ScrSize(4)/50 ...                        %  figure on the 
        ScrSize(3)/3.25 ScrSize(4)/1.10]);                      %  screen
    
        CPOSupTtl_Hndl = ...
            suptitle(['Patient Number ' num2str(Patient_Num(i))]);
        set(CPOSupTtl_Hndl,'FontSize',20,'FontWeight','bold')
        CPOLgnd_Strgs = cell(1,Num_RHCs(i));
        subplot(2,1,1)
        for j = 1:Num_RHCs(i)
            LV_CPOAve = mean(LV_CPO{i,j}(:));
            Days_PTx = Lngtdnl_OptimP{i}.data(1,j);
            plot(Days_PTx,LV_CPOAve, ...
                ['o' LnColor_Str{j}], ...
                'MarkerSize',9, ...
                'MarkerFaceColor',LnColor_Str{j}, ...
                'MarkerEdgeColor','k')
            hold on
            CPOLgnd_Strgs{j} = num2str(Days_PTx);
        end
        
        xlim([0 Lngtdnl_OptimP{i}.data(1,end)*1.15])
        if (ismember(Patient_Num(i),[363 558 572]))
            LgndLoc = 'SouthEast';
        else
            LgndLoc = 'NorthEast';
        end
        LegHndlCPOFigi_SP1 = legend(CPOLgnd_Strgs{:});          % Format subplot
        set(LegHndlCPOFigi_SP1,'Box','off', ...
            'FontSize',8,'Location',LgndLoc)
        set(gca,'FontSize',14, ...
            'FontWeight','bold','Box','off')
        xlabel('Days Post Transplant', ...
            'FontSize',18,'FontWeight','bold')
        ylabel('LV Power Output (W)', ...
            'FontSize',18,'FontWeight','bold')
        
        subplot(2,1,2)
        for j = 1:Num_RHCs(i)
            RV_CPOAve = mean(RV_CPO{i,j}(:));
            Days_PTx = Lngtdnl_OptimP{i}.data(1,j);
            plot(Days_PTx,RV_CPOAve, ...
                ['o' LnColor_Str{j}], ...
                'MarkerSize',9, ...
                'MarkerFaceColor',LnColor_Str{j}, ...
                'MarkerEdgeColor','k')
            hold on
        end
        
        xlim([0 Lngtdnl_OptimP{i}.data(1,end)*1.15])
        if (Patient_Num(i) == 0)
            LgndLoc = 'SouthEast';
        else
            LgndLoc = 'NorthEast';
        end
        if (Patient_Num(i) == 363)
            LgndLoc = 'SouthEast';
        else
            LgndLoc = 'NorthEast';
        end
        LegHndlCPOFigi_SP2 = legend(CPOLgnd_Strgs{:});          % Format subplot
        set(LegHndlCPOFigi_SP2,'Box','off', ...
            'FontSize',8,'Location',LgndLoc)
        set(gca,'FontSize',14, ...
            'FontWeight','bold','Box','off')
        xlabel('Days Post Transplant', ...
            'FontSize',18,'FontWeight','bold')
        ylabel('RV Power Output (W)', ...
            'FontSize',18,'FontWeight','bold')
        
    end
        
        
            
        