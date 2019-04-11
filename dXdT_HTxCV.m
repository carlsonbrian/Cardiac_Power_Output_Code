% ***********************************************************************************
%           d X d T   F U N C T I O N   for   H E A R T   T R A N S P L A N T
%               C A R D I O V A S C U L A R   S Y S T E M S   M O D E L
% ***********************************************************************************
%
%   This function contains the algebraic and differential expressions that describe
%   the reduced form cardiovascular systems model used in our patient specific
%   modeling of post heart transplant cardivascular function from right heart
%   catheterization measures
%
%   Model originally created on     17 January 2017
%   Model last modfied on           31 March   2019
%
%   Reproduced by       Payton Woodall, Amanda Colunga, 
%                       Mette Olufsen and Brian Carlson
%                       
%                       Physiological Systems Dynamics Laboratory
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
%
% ***********************************************************************************
%  START OF  	     d X d t  for   S M I T H   C A R D I O V A S C U L A R   
%                                         S Y S T E M S   M O D E L
% ***********************************************************************************

function [Var_Out] = dXdT_HTxCV(time,X,RHCParam_Struct,CVParam_Struct,varargin)

    % UNPACKING RHC PARAMETERS
    HR = RHCParam_Struct.HR;                    % Heart rate (beats/min)
    % UNPACKING MODEL PARAMETERS
    % Left ventricle free wall parameters
    E_es_lvf = CVParam_Struct.E_es_lvf;         % LV free wall elstnce (kPa/mL)
    V_d_lvf = CVParam_Struct.V_d_lvf;           % LV ES zero P volume (mL)
    P_0_lvf = CVParam_Struct.P_0_lvf;           % LV ED pressure param (kPa)
    lambda_lvf = CVParam_Struct.lambda_lvf;     % LV ED pressure param (1/mL)
    V_0_lvf = CVParam_Struct.V_0_lvf;           % LV ED pressure param (mL)
    % Right ventricle free wall parameters
    E_es_rvf = CVParam_Struct.E_es_rvf;         % RV free wall elstnce (kPa/mL)
    V_d_rvf = CVParam_Struct.V_d_rvf;           % RV ES zero P volume (mL)
    P_0_rvf = CVParam_Struct.P_0_rvf;           % RV ED pressure param (kPa)
    lambda_rvf = CVParam_Struct.lambda_rvf;     % RV ED pressure param (1/mL)
    V_0_rvf = CVParam_Struct.V_0_rvf;           % RV ED pressure param (mL)
    % Pulmonary artery and vein parameters
    E_es_pa = CVParam_Struct.E_es_pa;           % Pulm artery elastance (kPa/mL)
    E_es_pu = CVParam_Struct.E_es_pu;           % Pulm vein elastance (kPa/mL)
    R_pul = CVParam_Struct.R_pul;               % Pulm vasc resistnce (kPa*s/mL)
    P_th = CVParam_Struct.P_th;                 % Mean thoracic pressure (mmHg)
    % Aortic and vena cava parameters
    E_es_ao = CVParam_Struct.E_es_ao;           % Aorta elastance (kPa/mL)
    E_es_vc = CVParam_Struct.E_es_vc;           % Vena cava elastance (kPa/mL)
    R_sys = CVParam_Struct.R_sys;               % Syst art resistance (kPa*s/mL)
    % Heart valve parameters
    R_mt = CVParam_Struct.R_mt;                 % Mitral valve resist (kPa*s/mL)
    R_av = CVParam_Struct.R_av;                 % Aortic valve resist (kPa*s/mL)
    R_tc = CVParam_Struct.R_tc;                 % Tricspd vlv resist (kPa*s/mL)
    R_pv = CVParam_Struct.R_pv;                 % Pulmon vlv resist (kPa*s/mL)

    % UNPACK STATE VARIABLES
    V_lv = X(1);                                % Left ventricular volume (mL)
    V_rv = X(2);                                % Right ventricular volume (mL)
    V_pa = X(3);                                % Pulmonary artery volume (mL)
    V_pu = X(4);                                % Pulmonary vein volume (mL)
    V_ao = X(5);                                % Aortic/systemic volume (mL)
    V_vc = X(6);                                % Vena cava/venous volume (mL)
    
    % DRIVER FUNCTION
    T = 60 / HR;                                % Period (sec/beat)
    B = HR;                                     % Elastance function parameter
    C = T / 2;                                  % Elastance function parameter
    tau  = time - (floor(time/T) * T);
    e_t  = exp((-1) * B * (tau - C)^2);

    % LEFT VENTRICULAR PRESSURE
    P_es_lvf = E_es_lvf * (V_lv - V_d_lvf);
    P_ed_lvf = P_0_lvf * (exp(lambda_lvf * (V_lv - V_0_lvf)) - 1);
    P_lv = (e_t * P_es_lvf) + ((1-e_t) * P_ed_lvf) + P_th;
    
    % RIGHT VENTRICULAR PRESSURE
    P_es_rvf = E_es_rvf * (V_rv - V_d_rvf);
    P_ed_rvf = P_0_rvf * (exp(lambda_rvf * (V_rv - V_0_rvf)) - 1);
    P_rv = (e_t * P_es_rvf) + ((1-e_t) * P_ed_rvf) + P_th;
    
    % PULMONARY CIRCULATION PRESSURES AND FLOW
    P_pa = E_es_pa * (V_pa) + P_th;
    P_pu = E_es_pu * (V_pu) + P_th;
    Q_pul = (P_pa-P_pu) / R_pul;

    % SYSTEMIC CIRCULATION PRESSURES AND FLOW
    P_ao = E_es_ao * (V_ao);
    P_vc = E_es_vc * (V_vc);
    Q_sys = (P_ao-P_vc) / R_sys;
    
    % FLOWS THROUGH HEART VALVES
    if P_pu > P_lv
        Q_mt = (P_pu - P_lv)/R_mt;
    else
        Q_mt = 0;
    end
    if P_lv > P_ao
        Q_av = (P_lv - P_ao)/R_av;
    else
        Q_av = 0;
    end
    if P_vc > P_rv
        Q_tc = (P_vc - P_rv)/R_tc;
    else
        Q_tc = 0;
    end
    if P_rv > P_pa
        Q_pv = (P_rv - P_pa)/R_pv;
    else
        Q_pv = 0;
    end
    
    % DERIVATIVES OF STATE VARIABLES
    dVlvdt = Q_mt  - Q_av;
    dVrvdt = Q_tc  - Q_pv;
    dVpadt = Q_pv  - Q_pul;
    dVpudt = Q_pul - Q_mt;
    dVaodt = Q_av  - Q_sys;
    dVvcdt = Q_sys - Q_tc;

    % RATE OF CHANGE (RoC) OF STATE VARIABLES
    % Here the RoCs are returned unless we are running this 
    %  function outside the ODE solver to calculate other
    %  non-state variables to plot
    if (isempty(varargin))       
        RoC(1) = dVlvdt;
        RoC(2) = dVrvdt;
        RoC(3) = dVpadt;
        RoC(4) = dVpudt;
        RoC(5) = dVaodt;
        RoC(6) = dVvcdt;
        Var_Out = RoC';     
    else     
        CalcVars(1) = P_lv;
        CalcVars(2) = P_rv;
        CalcVars(3) = P_ao;
        CalcVars(4) = P_vc;
        CalcVars(5) = P_pa;
        CalcVars(6) = P_pu;
        Var_Out = CalcVars';
    end
    
end

