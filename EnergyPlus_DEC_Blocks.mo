within ;
package EnergyPlus_DEC_Blocks
  "Direct evaporative cooler with energyplus equation"

  package BaseClasses

    block Pad_Efficiency "To compute the evaporative cooling pad efficiency "
    //extends Buildings.BaseClasses.BaseIcon;
       Modelica.Blocks.Interfaces.RealInput vel_air "Velocity of the air (m/s)"
        annotation (Placement(transformation(extent={{-114,-48},{-74,-8}}),
            iconTransformation(extent={{-114,-48},{-74,-8}})));
      Modelica.Blocks.Interfaces.RealOutput efficiency "Efficiency of the cooling pad "
                                                       annotation (Placement(transformation(
              extent={{92,-12},{112,8}}), iconTransformation(extent={{92,-12},{112,8}})));
              parameter Real tk_pad = 1
                                       " thickness of the cooling pad in m";
              output Real Eff;
    equation
      Eff = 0.792714  + (0.958569*tk_pad) - (0.25193* vel_air) - (1.03215*(tk_pad^2)) + (0.0262659*(vel_air^2)) + (0.914869* (tk_pad*vel_air)) - (1.4821 *(vel_air*(tk_pad^2))) - 0.018992*((vel_air^3)*tk_pad^2) + (1.13137*((tk_pad^3)*vel_air)) + (0.0327622*(vel_air^3*tk_pad^2)) - (0.145384*(tk_pad^3)*(vel_air^2));
      efficiency = if Eff>1 then 1 elseif Eff<=0 then 0 else Eff;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Pad_Efficiency;

    block Water_consumption
       //extends Buildings.BaseClasses.BaseIcon;
      Modelica.Blocks.Interfaces.RealInput m_flo
        "mass flow rate of the water (kg/s)"
        annotation (Placement(transformation(extent={{-114,38},{-74,78}}),
            iconTransformation(extent={{-110,6},{-80,36}})));
      Modelica.Blocks.Interfaces.RealInput W_in
        "Humidity ratio of the inlet water (kg/kg)"
        annotation (Placement(transformation(extent={{-114,38},{-74,78}}),
            iconTransformation(extent={{-110,52},{-80,82}})));
      Modelica.Blocks.Interfaces.RealInput W_out
        "Humidity ratio of the outlet water(kg/kg)"
                                             annotation (Placement(transformation(
              extent={{-114,38},{-74,78}}), iconTransformation(extent={{-110,-40},{-80,
                -10}})));
      Modelica.Blocks.Interfaces.RealInput Density "Density of water (kg/m3)"
                                                                      annotation (
          Placement(transformation(extent={{-114,38},{-74,78}}), iconTransformation(
            extent={{-15,-15},{15,15}},
            rotation=0,
            origin={-95,-71})));


    Real W_drift "water consumption due to drift losses";
    Real W_Bl_dn " water consumption due to blowdown";
    parameter Real f_drift = 1 "Drift fraction";
    parameter Real Rcon= 1
                          "Ratio of solids in blowdown water";
                          //output Real Vol_water;

      Modelica.Blocks.Interfaces.RealOutput Vol_evap
        "evaporativel water consumption in m3/s" annotation (Placement(
            transformation(extent={{94,-6},{114,14}}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={100,8})));
      Modelica.Blocks.Interfaces.RealOutput Vol_water
        "Total volume of water consumed" annotation (Placement(transformation(
              extent={{88,-46},{108,-26}}),iconTransformation(extent={{90,-44},
                {114,-20}})));
    algorithm

      Vol_evap := (m_flo*abs(W_in - W_out))/Density;
      W_drift :=Vol_evap*f_drift;
      W_Bl_dn :=(Vol_evap/(Rcon - 1)) - W_drift;
      Vol_water :=Vol_evap + W_drift + W_Bl_dn;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Water_consumption;

    block OutletConditions
       //extends Buildings.BaseClasses.BaseIcon;
      Modelica.Blocks.Interfaces.RealInput Efficiency annotation (Placement(
            transformation(extent={{-116,-54},{-76,-14}}), iconTransformation(
              extent={{-116,-54},{-76,-14}})));
      Modelica.Blocks.Interfaces.RealInput Twb_air annotation (Placement(
            transformation(extent={{-114,32},{-74,72}}), iconTransformation(extent={
                {-116,-4},{-76,36}})));
      Modelica.Blocks.Interfaces.RealInput Tdb_air annotation (Placement(
            transformation(extent={{-114,32},{-74,72}}), iconTransformation(extent={
                {-116,46},{-76,86}})));

      Modelica.Blocks.Interfaces.RealOutput Wx_out "Humidity ratio of the exit air"
        annotation (Placement(transformation(extent={{94,26},{114,46}}), iconTransformation(extent={{94,-2},{114,18}})));

         output Real TDew "Dew point temperature,K";
         output Modelica.SIunits.Temperature Twb_out "Wet Bulb temperature,K";
         parameter Modelica.SIunits.Pressure p_atm = 101325
                                                       "atmospheric pressure(in pascals)";
         output Modelica.SIunits.Pressure p_wat "partial vapour pressure";
         output Real phi "Relative humidity";
         output Real Tdb_out "outlet temeprature";
         output Real Tdb_outC " Conversion to C from K";
         output Real Twb_outC "Coversion to C from K";
         output Real ed, ew,B,e;
    algorithm

      //computing Tdb of the leaving air
      Tdb_out :=Tdb_air- (Efficiency*(Tdb_air-Twb_air));

       //Considering Twb_air(in) = Twb_out (psychrometrics of evaporative cooling), and Q = 0
     Twb_out := Twb_air;

      // converting to celcius
      Tdb_outC := Modelica.SIunits.Conversions.to_degC(Tdb_out);
      Twb_outC := Modelica.SIunits.Conversions.to_degC(Twb_air);

     //computing RH(phi) using Tdb_out, Twb_out and atmospheric pressure

    ed := 6.108 * exp((17.27*Tdb_outC)/(237.3+Tdb_outC))
                                                       "Saturation Vapor Pressure at Dry Bulb (mb)";
    ew := 6.108 * exp((17.27*Twb_outC)/(237.3+Twb_outC))
                                                       "Saturation Vapor Pressure at wet bulb Bulb (mb)";
    e := ew - (0.00066 * (1+0.00115*Twb_outC)*(Tdb_outC-Twb_outC)*(p_atm/100));
    phi := 100*(e/ed);
    B := log(e/6.108)/17.27;
    TDew := (237.3*B)/(1-B);
    p_wat := (6.11* 10^((7.5*TDew)/(237.3+TDew)))*100;

      //p_wat:= Buildings.Utilities.Psychrometrics.Functions.pW_TDewPoi_amb(Modelica.SIunits.Conversions.from_degC(TDew));
      Wx_out:= Buildings.Utilities.Psychrometrics.Functions.X_pW(p_wat);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end OutletConditions;

    block NTU_effectiveness
      "Calculating the NTU and Effectiveness for the DEC model"
      Modelica.Blocks.Interfaces.RealInput m2_flow "Mass flow rate of air"
        annotation (Placement(transformation(extent={{-116,50},{-76,90}}),
            iconTransformation(extent={{-112,40},{-82,70}})));
      Modelica.Blocks.Interfaces.RealInput T_wb1
        "inlet wet bulb temperature of air" annotation (Placement(transformation(
              extent={{-116,50},{-76,90}}), iconTransformation(extent={{-112,6},{-82,
                36}})));
      Modelica.Blocks.Interfaces.RealInput T_db1
        "Inlet Drybulb temperature of air"
        annotation (Placement(transformation(extent={{-116,50},{-76,90}}),
            iconTransformation(extent={{-112,-28},{-82,2}})));
      Modelica.Blocks.Interfaces.RealInput Tw1 "Water inlet temperature"
        annotation (Placement(transformation(extent={{-116,50},{-76,90}}),
            iconTransformation(extent={{-112,-64},{-82,-34}})));
      Modelica.Blocks.Interfaces.RealOutput T_db2 "Outlet drybulb temperature"
        annotation (Placement(transformation(extent={{92,46},{112,66}}),
            iconTransformation(extent={{92,46},{112,66}})));
      Modelica.Blocks.Interfaces.RealOutput Tw2 "Outlet water temperature"
        annotation (Placement(transformation(extent={{92,46},{112,66}}),
            iconTransformation(extent={{92,18},{112,38}})));
      Modelica.Blocks.Interfaces.RealOutput X_w2 "Outlet humidity ratio"
        annotation (Placement(transformation(extent={{92,46},{112,66}}),
            iconTransformation(extent={{92,-8},{112,12}})));

      Real gama "constant depending on pad material and pad configu- ration like ξ";
      Real Nu " Nusselt number";
      Real Re " Reynolds number";
      Real Pr
             " Prandtl number";
      Real Mu " Dynamic viscosity";
      Real Cpa  " Specific heat of air";
      parameter Real p_atm " Atmospheric pressure";
      parameter Real k "Thermal conductivity";
      parameter Real v_a "air velocity";
      Real rho " air density";
      parameter Real Thickness "thickness of the cooling pad";
      parameter Real Area_s " area of the cooling pad";
      Real l
            "Characteristic length";
      Real NTU " Net transfer unit";
      parameter Real V "volume flow rate";
      Real phi " Intermediate variable - relative humidity";
      Real Mr " mass transfer coefficient";
      Real alfa, beta;
      Real T_wb2, TDew, p_wat,ed,ew,e,B,RH;

      Modelica.Blocks.Interfaces.RealOutput hm "mass transfer coefficient "
        annotation (Placement(transformation(extent={{92,46},{112,66}}),
            iconTransformation(extent={{92,-40},{112,-20}})));
      Modelica.Blocks.Interfaces.RealOutput hc
        "convective heat transfer coefficient" annotation (Placement(transformation(
              extent={{92,46},{112,66}}), iconTransformation(extent={{92,-66},{112,-46}})));
    algorithm

      // calculating Nusselts number to compute hc, as Nu gives a relationship between conductive and convective heat transfer
    l:= (Area_s / Thickness);
     gama := 0.1 * l^0.12;
     rho := p_atm/287.08 * T_wb1;
     Re := Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers.ReynoldsNumber(v_a, rho, Mu, l);
     Pr := Modelica.Fluid.Dissipation.Utilities.Functions.General.PrandtlNumber(Cpa, Mu, k);
     Nu:= gama * Re^0.8 * Pr^0.33;

     hc := Nu * k;
     //Lewis relationship between convective heat transfer coefficient and mass transfer coefficient (Le ~ 1 for water vapour and air)
     hm := hc/Cpa;
     //Net transfer unit calculation
     NTU := (hc * Area_s)/(m2_flow * Cpa);
     // calculating the efficiency of the cooling pad
     // calculated using alfa, beta and Mr

     alfa := (Tw1 - T_wb1)/(T_db1 - T_wb1);
     beta := 0.44 * (alfa ^ 0.932)* (NTU^0.089)* (Mr^(-0.116));
     Tw2 := Tw1 - beta*(T_db1-T_wb1);
     Mr := (m2_flow/V);
     phi := 1- exp(-1.07*(NTU^0.295))* alfa ^(-0.556) * Mr^(-0.051);
     //η= ϕη_1
     //η_1=1-exp⁡[-1.037 NTU]
     //η=(T_db1- T_db2)/(T_db1- T_wb1)

     T_db2 := T_db1 - ((T_db1 - T_wb1)*phi*(1-exp((-1.037)*NTU)));
     // Humidity ratio for the outlet condition
     //assuming WBT in = WBT out ;

     T_wb2:= T_wb1;

     //computing RH(phi) using Tdb_out, Twb_out and atmospheric pressure

    ed := 6.108 * exp((17.27*T_db2)/(237.3+T_db2))     "Saturation Vapor Pressure at Dry Bulb (mb)";
    ew := 6.108 * exp((17.27*T_wb2)/(237.3+T_wb2))     "Saturation Vapor Pressure at wet bulb Bulb (mb)";
    e := ew - (0.00066 * (1+0.00115*T_wb2)*(T_db2-T_wb2)*(p_atm/100));
    RH := 100*(e/ed);
    B := log(e/6.108)/17.27;
    TDew := (237.3*B)/(1-B);
    p_wat := (6.11* 10^((7.5*TDew)/(237.3+TDew)))*100;
    X_w2 := Buildings.Utilities.Psychrometrics.Functions.X_pW(p_wat);


      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end NTU_effectiveness;

    block HeatTransfer
      Modelica.Blocks.Interfaces.RealInput hc "Heat transfer coefficient"
        annotation (Placement(transformation(extent={{-140,80},{-100,120}}),
            iconTransformation(extent={{-132,100},{-102,130}})));
      Modelica.Blocks.Interfaces.RealInput hm "mass transfer coefficient"
        annotation (Placement(transformation(extent={{-140,-40},{-100,0}}),
            iconTransformation(extent={{-132,-36},{-100,-4}})));
      Modelica.Blocks.Interfaces.RealInput T_db1 "Inlet drybulb temperature"
        annotation (Placement(transformation(extent={{-140,40},{-100,80}}),
            iconTransformation(extent={{-132,64},{-100,96}})));
      Modelica.Blocks.Interfaces.RealInput T_db2 "Outlet drybulb temperature"
        annotation (Placement(transformation(extent={{-140,-80},{-100,-40}}),
            iconTransformation(extent={{-132,-72},{-100,-40}})));
      Modelica.Blocks.Interfaces.RealInput Wx_1 "Inlet humidity ratio" annotation (
          Placement(transformation(extent={{-140,0},{-100,40}}), iconTransformation(
              extent={{-132,28},{-100,60}})));
      Modelica.Blocks.Interfaces.RealInput Wx_2 "Outlet humidity ratio" annotation (
         Placement(transformation(extent={{-140,-118},{-100,-78}}),
                                                                iconTransformation(
              extent={{-132,-110},{-100,-78}})));
      Modelica.Blocks.Interfaces.RealOutput Q_total "Total heat transfer"
        annotation (Placement(transformation(extent={{120,-10},{140,10}}),
            iconTransformation(extent={{120,-10},{140,10}})));
      Modelica.Blocks.Interfaces.RealOutput Q_latent "Latent heat transfer"
        annotation (Placement(transformation(extent={{120,-70},{140,-50}}),
            iconTransformation(extent={{120,-70},{140,-50}})));

            parameter Real Area_s= 1 "Area of the cooling pad";
            parameter Real hvs= 2675.43*1000 "specific enthalphy of vapour at surface temperature";
      Modelica.Blocks.Interfaces.RealOutput Q_sensible
                                                     "sensible heat transfer"
        annotation (Placement(transformation(extent={{120,50},{140,70}}),
            iconTransformation(extent={{120,50},{140,70}})));
    algorithm
      Q_sensible := hc * Area_s *(T_db2- T_db1);
      Q_latent := hvs * hm * Area_s *( Wx_2 - Wx_1);
      Q_total := Q_sensible + Q_latent;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-120},
                {120,120}})),                                        Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-120},{120,120}})));
    end HeatTransfer;

    block WaterConsumption_Detailed
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end WaterConsumption_Detailed;

    model NTU_Eff_calculation
      Modelica.Blocks.Interfaces.RealInput ma_flow "Mass flow rate of air,kg/s"
        annotation (Placement(transformation(extent={{-116,50},{-76,90}}),
            iconTransformation(extent={{-112,40},{-82,70}})));
      Modelica.Blocks.Interfaces.RealInput Ta_wb1
        "inlet wet bulb temperature of air,C" annotation (Placement(transformation(
              extent={{-116,50},{-76,90}}), iconTransformation(extent={{-112,6},{-82,
                36}})));
      Modelica.Blocks.Interfaces.RealInput Ta_db1
        "Inlet Drybulb temperature of air,C" annotation (Placement(transformation(
              extent={{-116,50},{-76,90}}), iconTransformation(extent={{-112,-28},{-82,
                2}})));
      Modelica.Blocks.Interfaces.RealInput Tw1 "Water inlet temperature,C"
        annotation (Placement(transformation(extent={{-116,50},{-76,90}}),
            iconTransformation(extent={{-112,-64},{-82,-34}})));
      Modelica.Blocks.Interfaces.RealOutput Ta_db2 "Outlet drybulb temperature,C"
        annotation (Placement(transformation(extent={{92,32},{112,52}}),
            iconTransformation(extent={{92,32},{112,52}})));
      Modelica.Blocks.Interfaces.RealOutput X_2 "Outlet humidity ratio , kg/kg" annotation (
         Placement(transformation(extent={{92,58},{112,78}}), iconTransformation(
              extent={{92,58},{112,78}})));
      Modelica.Blocks.Interfaces.RealOutput hm "mass transfer coefficient "
        annotation (Placement(transformation(extent={{92,46},{112,66}}),
            iconTransformation(extent={{92,6},{112,26}})));
      Modelica.Blocks.Interfaces.RealOutput hc
        "convective heat transfer coefficient" annotation (Placement(transformation(
              extent={{92,46},{112,66}}), iconTransformation(extent={{92,-66},{112,-46}})));
      Modelica.Blocks.Interfaces.RealInput v_a
        "Velocity of air passing through the cooling pad , m/s" annotation (Placement(
            transformation(extent={{-116,50},{-76,90}}), iconTransformation(
            extent={{-15,15},{15,-15}},
            rotation=-90,
            origin={-35,101})));
      Modelica.Blocks.Interfaces.RealInput mw_flow "Massflow rate of water, kg/s"
        annotation (Placement(transformation(extent={{-116,50},{-76,90}}),
            iconTransformation(
            extent={{-15,15},{15,-15}},
            rotation=90,
            origin={-31,-95})));

    //parameters

      parameter Real p_atm = 101325 " Atmospheric pressure, pa ";
      parameter Real k = 0.045 "Thermal conductivity";
      parameter Real Thickness= 0.15 "Thickness of the cooling pad,m";
      parameter Real length = 2.1 "Length of the cooling pad,m";
      parameter Real height = 0.91 "Height of the cooling pad,m";
      parameter Real p_V = 440 "specific contact area, m2/m3";

     //Variables

      Real Area_s= length * height
                                  " area of the cooling pad, m2";
      Real V= Area_s*Thickness "volume of the cooling pad,m3";
      Real TSA = (2*Area_s)+(2*Thickness*(length+height))
                                                         "total surface area of the cooling pad,m2";
      Real Nu " Nusselt number";
      Real Re " Reynolds number";
      Real Pr " Prandtl number";
      Real Mu = 8* (10^(-6))*exp(0.0027*Ta_db1)  " Dynamic viscosity, mpa s";
      Real Cpa = 1000 " Specific heat of air J/kg K";
      Real rho " air density, kg/m3";
      Real l= (V/ TSA)
                      "Characteristic length,m";
      Real NTU " Net transfer unit";
      Real eff,B,Ta_db2C, Ta_wb2C "intermediate variables";
      Real Ta_wb2 "Wet bulb temperature of the outlet air";
      Real TDew "Dew point temperature of the outlet air";
      Real p_wat "Partial pressure";
      Real ed "Saturation Vapor Pressure at Dry Bulb (mb)";
      Real ew "Saturation Vapor Pressure at Wet Bulb (mb)";
      Real e "Actual Vapor Pressure (mb)";
      Real RH "Relative humidity, %";
      Real ka = 0.02706 "Thermal conductivity of air, W/mK";


      Modelica.Blocks.Interfaces.RealOutput DeltaP "Pressure drop (pascal)"
        annotation (Placement(transformation(extent={{92,46},{112,66}}),
            iconTransformation(extent={{92,-92},{112,-72}})));


    equation

      //calculating Nusselts number to compute hc,hm and NTU

      rho =  p_atm/(287.08 * (Ta_wb1));
      Re =  Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers.ReynoldsNumber(v_a, rho, Mu, l);
      Pr =  Modelica.Fluid.Dissipation.Utilities.Functions.General.PrandtlNumber(Cpa, Mu, ka);
      Nu=  0.192*((l/Thickness)^(-0.191))*(((Ta_db1- Ta_wb1)/Tw1)^(-0.03))* (Re^0.682) * (Pr^0.33);
      hc =  (Nu * k)/l;
      NTU =  (hc * TSA)/(ma_flow * Cpa);

     //Lewis relationship between convective heat transfer coefficient and mass transfer coefficient (Le ~ 1 for water vapour and air)
      hm =  hc/Cpa;
      eff =  1- exp((-hc*p_V*Thickness)/(v_a*rho*Cpa));
      Ta_db2 = Ta_db1 - eff*(Ta_db1 - Ta_wb1);

     //assuming WBT in = WBT out
      Ta_wb2=  Ta_wb1;

      Ta_db2C =   Modelica.SIunits.Conversions.to_degC(Ta_db2);
      Ta_wb2C =  Modelica.SIunits.Conversions.to_degC(Ta_wb2);

     //computing RH and Humidity ratio X_2 using Tdb_out, Twb_out and atmospheric pressure

     ed =  6.108 * exp((17.27*Ta_db2C)/(237.3+Ta_db2C));
     ew =  6.108 * exp((17.27*Ta_wb2C)/(237.3+Ta_wb2C));
     e =  ew - (0.00066 * (1+0.00115*Ta_wb2C)*(Ta_db2C-Ta_wb2C)*(p_atm/100));
     RH =  100*(e/ed);
     B =  log(e/6.108)/17.27;
     TDew =  (237.3*B)/(1-B);
     p_wat =  (6.11* 10^((7.5*(TDew))/(237.3+(TDew))))*100;
     X_2 =  Buildings.Utilities.Psychrometrics.Functions.X_pW(p_wat);


    // pressure drop module
    DeltaP = 0.124*((l/Thickness)^(-1.033))*(1 + 1825*(mw_flow/ma_flow)*((rho*v_a^2)
        /2));


       annotation (Line(points={{102,56},{58,56},{58,56},{102,56}},
            color={0,0,127}),
                  Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end NTU_Eff_calculation;

    record Breezair_icon_210_RPM560 "System performance curves for Breezair icon 210 under speed 9 (RPM 560)"
      extends Buildings.Fluid.Movers.Data.Generic(
        use_powerCharacteristic=true,
        speed_rpm_nominal=560,
        power(V_flow={1.53936,1.60388,1.68131,1.77380,1.88134,1.98242,2.07276,2.12223}, P={1059.89,1069.94,1090.034,1112.14,1140.27,1162.39,1186.50,1200.56}),
        pressure(V_flow={1.534,1.597,1.653,1.730,1.808,1.874,1.944,2.021,2.092,2.126}, dp={174.188,159.2576,144.3272,121.9316,99.536,79.6288,57.2332,34.8376,12.422}));
    end Breezair_icon_210_RPM560;

    model FourPortHeatMassExchanger
      "Model transporting two fluid streams between four ports with storing mass or energy"
      extends Buildings.Fluid.Interfaces.PartialFourPortInterface(
        port_a1(h_outflow(start=h1_outflow_start)),
        port_b1(h_outflow(start=h1_outflow_start)),
        port_a2(h_outflow(start=h2_outflow_start)),
        port_b2(h_outflow(start=h2_outflow_start)));
      extends Buildings.Fluid.Interfaces.FourPortFlowResistanceParameters(
         final computeFlowResistance1=true, final computeFlowResistance2=true);

      parameter Modelica.SIunits.Time tau1 = 30 "Time constant at nominal flow"
         annotation (Dialog(tab = "Dynamics", group="Nominal condition"));
      parameter Modelica.SIunits.Time tau2 = 30 "Time constant at nominal flow"
         annotation (Dialog(tab = "Dynamics", group="Nominal condition"));

      // Advanced
      parameter Boolean homotopyInitialization = true "= true, use homotopy method"
        annotation(Evaluate=true, Dialog(tab="Advanced"));

      // Assumptions
      parameter Modelica.Fluid.Types.Dynamics energyDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial
        "Type of energy balance: dynamic (3 initialization options) or steady state"
        annotation(Evaluate=true, Dialog(tab = "Dynamics", group="Equations"));
      parameter Modelica.Fluid.Types.Dynamics massDynamics=energyDynamics
        "Type of mass balance: dynamic (3 initialization options) or steady state"
        annotation(Evaluate=true, Dialog(tab = "Dynamics", group="Equations"));

      // Initialization
      parameter Medium1.AbsolutePressure p1_start = Medium1.p_default
        "Start value of pressure"
        annotation(Dialog(tab = "Initialization", group = "Medium 1"));
      parameter Medium1.Temperature T1_start = Medium1.T_default
        "Start value of temperature"
        annotation(Dialog(tab = "Initialization", group = "Medium 1"));
      parameter Medium1.MassFraction X1_start[Medium1.nX] = Medium1.X_default
        "Start value of mass fractions m_i/m"
        annotation (Dialog(tab="Initialization", group = "Medium 1", enable=Medium1.nXi > 0));
      parameter Medium1.ExtraProperty C1_start[Medium1.nC](
        final quantity=Medium1.extraPropertiesNames)=fill(0, Medium1.nC)
        "Start value of trace substances"
        annotation (Dialog(tab="Initialization", group = "Medium 1", enable=Medium1.nC > 0));
      parameter Medium1.ExtraProperty C1_nominal[Medium1.nC](
        final quantity=Medium1.extraPropertiesNames) = fill(1E-2, Medium1.nC)
        "Nominal value of trace substances. (Set to typical order of magnitude.)"
       annotation (Dialog(tab="Initialization", group = "Medium 1", enable=Medium1.nC > 0));

      parameter Medium2.AbsolutePressure p2_start = Medium2.p_default
        "Start value of pressure"
        annotation(Dialog(tab = "Initialization", group = "Medium 2"));
      parameter Medium2.Temperature T2_start = Medium2.T_default
        "Start value of temperature"
        annotation(Dialog(tab = "Initialization", group = "Medium 2"));
      parameter Medium2.MassFraction X2_start[Medium2.nX] = Medium2.X_default
        "Start value of mass fractions m_i/m"
        annotation (Dialog(tab="Initialization", group = "Medium 2", enable=Medium2.nXi > 0));
      parameter Medium2.ExtraProperty C2_start[Medium2.nC](
        final quantity=Medium2.extraPropertiesNames)=fill(0, Medium2.nC)
        "Start value of trace substances"
        annotation (Dialog(tab="Initialization", group = "Medium 2", enable=Medium2.nC > 0));
      parameter Medium2.ExtraProperty C2_nominal[Medium2.nC](
        final quantity=Medium2.extraPropertiesNames) = fill(1E-2, Medium2.nC)
        "Nominal value of trace substances. (Set to typical order of magnitude.)"
       annotation (Dialog(tab="Initialization", group = "Medium 2", enable=Medium2.nC > 0));

      Modelica.SIunits.HeatFlowRate Q1_flow = vol1.heatPort.Q_flow
        "Heat flow rate into medium 1";
      Modelica.SIunits.HeatFlowRate Q2_flow = vol2.heatPort.Q_flow
        "Heat flow rate into medium 2";

      replaceable Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatPort vol1
        constrainedby
        Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatPort(
            redeclare final package Medium = Medium1,
            nPorts = 2,
            V=m1_flow_nominal*tau1/rho1_nominal,
            final allowFlowReversal=allowFlowReversal1,
            final m_flow_nominal=m1_flow_nominal,
            energyDynamics=if tau1 > Modelica.Constants.eps
                             then energyDynamics else
                             Modelica.Fluid.Types.Dynamics.SteadyState,
            massDynamics=if tau1 > Modelica.Constants.eps
                             then massDynamics else
                             Modelica.Fluid.Types.Dynamics.SteadyState,
            final p_start=p1_start,
            final T_start=T1_start,
            final X_start=X1_start,
            final C_start=C1_start,
            final C_nominal=C1_nominal,
            mSenFac=1) "Volume for fluid 1"
        annotation (Placement(transformation(extent={{-10,70}, {10,50}})));

      replaceable Buildings.Fluid.MixingVolumes.MixingVolume vol2
        constrainedby
        Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatPort(
            redeclare final package Medium = Medium2,
            nPorts = 2,
            V=m2_flow_nominal*tau2/rho2_nominal,
            final allowFlowReversal=allowFlowReversal2,
            mSenFac=1,
            final m_flow_nominal = m2_flow_nominal,
            energyDynamics=if tau2 > Modelica.Constants.eps
                             then energyDynamics else
                             Modelica.Fluid.Types.Dynamics.SteadyState,
            massDynamics=if tau2 > Modelica.Constants.eps
                             then massDynamics else
                             Modelica.Fluid.Types.Dynamics.SteadyState,
            final p_start=p2_start,
            final T_start=T2_start,
            final X_start=X2_start,
            final C_start=C2_start,
            final C_nominal=C2_nominal) "Volume for fluid 2"
       annotation (Placement(transformation(
            origin={-40,-60},
            extent={{-10,10},{10,-10}},
            rotation=180)));

      Buildings.Fluid.FixedResistances.PressureDrop preDro1(
        redeclare final package Medium = Medium1,
        final m_flow_nominal=m1_flow_nominal,
        final deltaM=deltaM1,
        final allowFlowReversal=allowFlowReversal1,
        final show_T=false,
        final from_dp=from_dp1,
        final linearized=linearizeFlowResistance1,
        final homotopyInitialization=homotopyInitialization,
        final dp_nominal=dp1_nominal) "Flow resistance of fluid 1"
        annotation (Placement(transformation(extent={{-60,70},{-40,90}})));

      Buildings.Fluid.FixedResistances.PressureDrop preDro2(
        redeclare final package Medium = Medium2,
        final m_flow_nominal=m2_flow_nominal,
        final deltaM=deltaM2,
        final allowFlowReversal=allowFlowReversal2,
        final show_T=false,
        final from_dp=from_dp2,
        final linearized=linearizeFlowResistance2,
        final homotopyInitialization=homotopyInitialization,
        final dp_nominal=dp2_nominal) "Flow resistance of fluid 2"
        annotation (Placement(transformation(extent={{0,-90},{-20,-70}})));

      Buildings.Fluid.Sensors.Temperature senTem_b2(redeclare package Medium =
            Medium2)
        annotation (Placement(transformation(extent={{-96,-60},{-76,-40}})));
      Buildings.Fluid.Sensors.Temperature senTem_a2(redeclare package Medium =
            Medium2)
        annotation (Placement(transformation(extent={{70,-60},{90,-40}})));
      Buildings.Fluid.Sensors.Temperature senTem_b1(redeclare package Medium =
            Medium1)
        annotation (Placement(transformation(extent={{70,60},{90,80}})));
      Buildings.Fluid.Sensors.Temperature senTem_a1(redeclare package Medium =
            Medium1)
        annotation (Placement(transformation(extent={{-98,80},{-78,100}})));
      Buildings.Fluid.Sensors.MassFlowRate senMasFlo1(redeclare package Medium
          = Medium1)
        annotation (Placement(transformation(extent={{-84,70},{-64,90}})));
      Buildings.Fluid.Sensors.MassFlowRate senMasFlo2(redeclare package Medium
          = Medium2)
        annotation (Placement(transformation(extent={{90,-90},{70,-70}})));
      Buildings.Fluid.Sensors.Velocity senVel(
        redeclare package Medium = Medium2,
        m_flow_nominal=m2_flow_nominal,
        A=A) annotation (Placement(transformation(extent={{64,-90},{44,-70}})));
      parameter Modelica.SIunits.Area A "Cross sectional area of flow channel";
      Buildings.Fluid.Sensors.TemperatureWetBulbTwoPort senWetBul(redeclare
          package Medium = Medium2, m_flow_nominal=m2_flow_nominal)
        annotation (Placement(transformation(extent={{38,-90},{18,-70}})));
    protected
      parameter Medium1.ThermodynamicState sta1_nominal=Medium1.setState_pTX(
          T=Medium1.T_default, p=Medium1.p_default, X=Medium1.X_default);
      parameter Modelica.SIunits.Density rho1_nominal=Medium1.density(sta1_nominal)
        "Density, used to compute fluid volume";
      parameter Medium2.ThermodynamicState sta2_nominal=Medium2.setState_pTX(
          T=Medium2.T_default, p=Medium2.p_default, X=Medium2.X_default);
      parameter Modelica.SIunits.Density rho2_nominal=Medium2.density(sta2_nominal)
        "Density, used to compute fluid volume";

      parameter Medium1.ThermodynamicState sta1_start=Medium1.setState_pTX(
          T=T1_start, p=p1_start, X=X1_start);
      parameter Modelica.SIunits.SpecificEnthalpy h1_outflow_start = Medium1.specificEnthalpy(sta1_start)
        "Start value for outflowing enthalpy";
      parameter Medium2.ThermodynamicState sta2_start=Medium2.setState_pTX(
          T=T2_start, p=p2_start, X=X2_start);
      parameter Modelica.SIunits.SpecificEnthalpy h2_outflow_start = Medium2.specificEnthalpy(sta2_start)
        "Start value for outflowing enthalpy";

    initial equation
      // Check for tau1
      assert((energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) or
              tau1 > Modelica.Constants.eps,
    "The parameter tau1, or the volume of the model from which tau may be derived, is unreasonably small.
 You need to set energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState to model steady-state.
 Received tau1 = "     + String(tau1) + "\n");
      assert((massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) or
              tau1 > Modelica.Constants.eps,
    "The parameter tau1, or the volume of the model from which tau may be derived, is unreasonably small.
 You need to set massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState to model steady-state.
 Received tau1 = "     + String(tau1) + "\n");

     // Check for tau2
      assert((energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) or
              tau2 > Modelica.Constants.eps,
    "The parameter tau2, or the volume of the model from which tau may be derived, is unreasonably small.
 You need to set energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState to model steady-state.
 Received tau2 = "     + String(tau2) + "\n");
      assert((massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) or
              tau2 > Modelica.Constants.eps,
    "The parameter tau2, or the volume of the model from which tau may be derived, is unreasonably small.
 You need to set massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState to model steady-state.
 Received tau2 = "     + String(tau2) + "\n");

    equation
      connect(vol1.ports[2], port_b1) annotation (Line(
          points={{0,70},{20,70},{20,60},{100,60}},
          color={0,127,255}));
      connect(vol2.ports[2], port_b2) annotation (Line(
          points={{-40,-70},{-60,-70},{-60,-60},{-100,-60}},
          color={0,127,255}));
      connect(preDro1.port_b, vol1.ports[1]) annotation (Line(
          points={{-40,80},{0,80},{0,70}},
          color={0,127,255}));
      connect(preDro2.port_b, vol2.ports[1]) annotation (Line(
          points={{-20,-80},{-40,-80},{-40,-70}},
          color={0,127,255}));
      connect(port_b2, senTem_b2.port)
        annotation (Line(points={{-100,-60},{-86,-60}}, color={0,127,255}));
      connect(port_a2, senTem_a2.port)
        annotation (Line(points={{100,-60},{80,-60}}, color={0,127,255}));
      connect(port_b1, senTem_b1.port)
        annotation (Line(points={{100,60},{80,60}}, color={0,127,255}));
      connect(port_a1, senTem_a1.port) annotation (Line(points={{-100,60},{-90,
              60},{-90,80},{-88,80}}, color={0,127,255}));
      connect(senTem_a1.port, senMasFlo1.port_a)
        annotation (Line(points={{-88,80},{-84,80}}, color={0,127,255}));
      connect(preDro1.port_a, senMasFlo1.port_b)
        annotation (Line(points={{-60,80},{-64,80}}, color={0,127,255}));
      connect(senMasFlo2.port_a, port_a2) annotation (Line(points={{90,-80},{92,
              -80},{92,-60},{100,-60}}, color={0,127,255}));
      connect(senMasFlo2.port_b, senVel.port_a)
        annotation (Line(points={{70,-80},{64,-80}}, color={0,127,255}));
      connect(senWetBul.port_a, senVel.port_b)
        annotation (Line(points={{38,-80},{44,-80}}, color={0,127,255}));
      connect(senWetBul.port_b, preDro2.port_a)
        annotation (Line(points={{18,-80},{0,-80}}, color={0,127,255}));
      annotation (
        Documentation(info="<html>
<p>
This component transports two fluid streams between four ports.
It provides the basic model for implementing a dynamic heat exchanger.
</p>
<p>
The model can be used as-is, although there will be no heat or mass transfer
between the two fluid streams.
To add heat transfer, heat flow can be added to the heat port of the two volumes.
See for example
<a href=\"Buildings.Fluid.Chillers.Carnot_y\">
Buildings.Fluid.Chillers.Carnot_y</a>.
To add moisture input into (or moisture output from) volume <code>vol2</code>,
the model can be replaced with
<a href=\"modelica://Buildings.Fluid.MixingVolumes.MixingVolumeMoistAir\">
Buildings.Fluid.MixingVolumes.MixingVolumeMoistAir</a>.
</p>
<h4>Implementation</h4>
<p>
The variable names follow the conventions used in
<a href=\"modelica://Modelica.Fluid.Examples.HeatExchanger.BaseClasses.BasicHX\">
Modelica.Fluid.Examples.HeatExchanger.BaseClasses.BasicHX</a>.
</p>
</html>",     revisions="<html>
<ul>
<li>
October 23, 2017, by Michael Wetter:<br/>
Made volume <code>vol1</code> replaceable. This is required for
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/1013\">Buildings, issue 1013</a>.
</li>
<li>
December 1, 2016, by Michael Wetter:<br/>
Updated model as <code>use_dh</code> is no longer a parameter in the pressure drop model.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/480\">#480</a>.
</li>
<li>
April 11, 2016 by Michael Wetter:<br/>
Corrected wrong hyperlink in documentation for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/450\">issue 450</a>.
</li>
<li>
January 26, 2016, by Michael Wetter:<br/>
Set <code>quantity</code> attributes.
</li>
<li>
November 13, 2015, by Michael Wetter:<br/>
Changed assignments of start values in <code>extends</code> statement.
This is for issue
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/299\">#299</a>.
</li>
<li>
June 2, 2015, by Filip Jorissen:<br/>
Removed final modifier from <code>mSenFac</code> in
<code>vol1</code> and <code>vol2</code>.
This is for issue
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/258=\">#258</a>.
</li>
<li>
May 6, 2015, by Michael Wetter:<br/>
Added missing propagation of <code>allowFlowReversal</code> to
instances <code>vol1</code> and <code>vol2</code>.
This is for issue
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/412\">#412</a>.
</li>
<li>
October 6, 2014, by Michael Wetter:<br/>
Changed medium declaration in pressure drop elements to be final.
</li>
<li>
May 28, 2014, by Michael Wetter:<br/>
Removed <code>annotation(Evaluate=true)</code> for parameters <code>tau1</code>
and <code>tau2</code>.
This is needed to allow changing the time constant after translation.
</li>
<li>
November 12, 2013, by Michael Wetter:<br/>
Removed <code>import Modelica.Constants</code> statement.
</li>
<li>
October 8, 2013, by Michael Wetter:<br/>
Removed parameter <code>show_V_flow</code>.
</li>
<li>
September 26, 2013, by Michael Wetter:<br/>
Removed unrequired <code>sum</code> operator.
</li>
<li>
February 6, 2012, by Michael Wetter:<br/>
Updated documentation.
</li>
<li>
February 3, 2012, by Michael Wetter:<br/>
Removed assignment of <code>m_flow_small</code> as it is no
longer used in its base class.
</li>
<li>
July 29, 2011, by Michael Wetter:
<ul>
<li>
Changed values of
<code>h_outflow_a1_start</code>,
<code>h_outflow_b1_start</code>,
<code>h_outflow_a2_start</code> and
<code>h_outflow_b2_start</code>, and
declared them as final.
</li>
<li>
Set nominal values for <code>vol1.C</code> and <code>vol2.C</code>.
</li>
</ul>
</li>
<li>
July 11, 2011, by Michael Wetter:<br/>
Changed parameterization of fluid volume so that steady-state balance is
used when <code>tau = 0</code>.
</li>
<li>
March 25, 2011, by Michael Wetter:<br/>
Added homotopy operator.
</li>
<li>
April 13, 2009, by Michael Wetter:<br/>
Added model to compute flow friction.
</li>
<li>
September 10, 2008 by Michael Wetter:<br/>
Added <code>stateSelect=StateSelect.always</code> for temperature of volume 1.
</li>
<li>
Changed temperature sensor from Celsius to Kelvin.
Unit conversion should be made during output
processing.
<li>
August 5, 2008, by Michael Wetter:<br/>
Replaced instances of <code>Delays.DelayFirstOrder</code> with instances of
<code>MixingVolumes.MixingVolume</code>. This allows to extract liquid for a condensing cooling
coil model.
</li>
<li>
March 25, 2008, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"),
        Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics={
            Rectangle(
              extent={{-70,80},{70,-80}},
              lineColor={0,0,255},
              pattern=LinePattern.None,
              fillColor={95,95,95},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-99,64},{102,54}},
              lineColor={0,0,255},
              pattern=LinePattern.None,
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-99,-56},{102,-66}},
              lineColor={0,0,255},
              pattern=LinePattern.None,
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid)}));
    end FourPortHeatMassExchanger;
  end BaseClasses;

  model DEC_EP_Method1

        extends Buildings.Fluid.Interfaces.PartialFourPortInterface(
        replaceable package Medium2 =
          Modelica.Media.Interfaces.PartialCondensingGases,
      allowFlowReversal2=false,
      allowFlowReversal1=false,
      port_a1(h_outflow(start=h1_outflow_start)),
      port_b1(h_outflow(start=h1_outflow_start)),
      port_a2(h_outflow(start=h2_outflow_start)),
      port_b2(h_outflow(start=h2_outflow_start)));
    extends Buildings.Fluid.Interfaces.FourPortFlowResistanceParameters(
       final computeFlowResistance1=true, final computeFlowResistance2=true);

  parameter Modelica.SIunits.Thickness PadThickness=0.3
      "Evaporative cooling pad thickness(m)"
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

  parameter Modelica.SIunits.Area PadArea= 0.2
      "Evaporative cooling pad area(m2)"
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

   parameter Real DriftFactor = 0.2
      "Drift factor[user input]"
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

   parameter Real Rcon = 0.2
      "Ratio of solids in the blowdown water[user input]"
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

    parameter Modelica.SIunits.Time tau1 = 1 "Time constant at nominal flow"
       annotation (Dialog(tab = "Dynamics", group="Nominal condition"));
    parameter Modelica.SIunits.Time tau2 = 1 "Time constant at nominal flow"
       annotation (Dialog(tab = "Dynamics", group="Nominal condition"));

    // Advanced
    parameter Boolean homotopyInitialization = true "= true, use homotopy method"
      annotation(Evaluate=true, Dialog(tab="Advanced"));

    // Assumptions
    parameter Modelica.Fluid.Types.Dynamics energyDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial
      "Type of energy balance: dynamic (3 initialization options) or steady state"
      annotation(Evaluate=true, Dialog(tab = "Dynamics", group="Equations"));
    parameter Modelica.Fluid.Types.Dynamics massDynamics=if tau1 > Modelica.Constants.eps
         then energyDynamics else Modelica.Fluid.Types.Dynamics.SteadyState
      "Type of mass balance: dynamic (3 initialization options) or steady state"
      annotation(Evaluate=true, Dialog(tab = "Dynamics", group="Equations"));

    // Initialization
    parameter Medium1.AbsolutePressure p1_start = Medium1.p_default
      "Start value of pressure"
      annotation(Dialog(tab = "Initialization", group = "Medium 1"));
    parameter Medium1.Temperature T1_start = Medium1.T_default
      "Start value of temperature"
      annotation(Dialog(tab = "Initialization", group = "Medium 1"));
    parameter Medium1.MassFraction X1_start[Medium1.nX] = Medium1.X_default
      "Start value of mass fractions m_i/m"
      annotation (Dialog(tab="Initialization", group = "Medium 1", enable=Medium1.nXi > 0));
    parameter Medium1.ExtraProperty C1_start[Medium1.nC](
      final quantity=Medium1.extraPropertiesNames)=fill(0, Medium1.nC)
      "Start value of trace substances"
      annotation (Dialog(tab="Initialization", group = "Medium 1", enable=Medium1.nC > 0));
    parameter Medium1.ExtraProperty C1_nominal[Medium1.nC](
      final quantity=Medium1.extraPropertiesNames) = fill(1E-2, Medium1.nC)
      "Nominal value of trace substances. (Set to typical order of magnitude.)"
     annotation (Dialog(tab="Initialization", group = "Medium 1", enable=Medium1.nC > 0));

    parameter Medium2.AbsolutePressure p2_start = Medium2.p_default
      "Start value of pressure"
      annotation(Dialog(tab = "Initialization", group = "Medium 2"));
    parameter Medium2.Temperature T2_start = Medium2.T_default
      "Start value of temperature"
      annotation(Dialog(tab = "Initialization", group = "Medium 2"));
    parameter Medium2.MassFraction X2_start[Medium2.nX] = Medium2.X_default
      "Start value of mass fractions m_i/m"
      annotation (Dialog(tab="Initialization", group = "Medium 2", enable=Medium2.nXi > 0));
    parameter Medium2.ExtraProperty C2_start[Medium2.nC](
      final quantity=Medium2.extraPropertiesNames)=fill(0, Medium2.nC)
      "Start value of trace substances"
      annotation (Dialog(tab="Initialization", group = "Medium 2", enable=Medium2.nC > 0));
    parameter Medium2.ExtraProperty C2_nominal[Medium2.nC](
      final quantity=Medium2.extraPropertiesNames) = fill(1E-2, Medium2.nC)
      "Nominal value of trace substances. (Set to typical order of magnitude.)"
     annotation (Dialog(tab="Initialization", group = "Medium 2", enable=Medium2.nC > 0));

    Modelica.SIunits.HeatFlowRate Q1_flow = vol1.heatPort.Q_flow
      "Heat flow rate into medium 1";
    Modelica.SIunits.HeatFlowRate Q2_flow = vol2.heatPort.Q_flow
      "Heat flow rate into medium 2";

    replaceable Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatPort vol1(nPorts=2)
      constrainedby
      Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatPort(
          redeclare final package Medium = Medium1,
          nPorts = 2,
          V=m1_flow_nominal*tau1/rho1_nominal,
          final allowFlowReversal=allowFlowReversal1,
          final m_flow_nominal=m1_flow_nominal,
          energyDynamics=if tau1 > Modelica.Constants.eps
                           then energyDynamics else
                           Modelica.Fluid.Types.Dynamics.SteadyState,
          massDynamics=if tau1 > Modelica.Constants.eps
                           then massDynamics else
                           Modelica.Fluid.Types.Dynamics.SteadyState,
          final p_start=p1_start,
          final T_start=T1_start,
          final X_start=X1_start,
          final C_start=C1_start,
          final C_nominal=C1_nominal,
          mSenFac=1) "Volume for fluid 1"
      annotation (Placement(transformation(extent={{26,60},{46,80}})));

    replaceable Buildings.Fluid.MixingVolumes.MixingVolumeMoistAir
                                                           vol2(nPorts=2)
      constrainedby
      Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatPort(
          redeclare final package Medium = Medium2,
          nPorts = 2,
          V=m2_flow_nominal*tau2/rho2_nominal,
          final allowFlowReversal=allowFlowReversal2,
          mSenFac=1,
          final m_flow_nominal = m2_flow_nominal,
          energyDynamics=if tau2 > Modelica.Constants.eps
                           then energyDynamics else
                           Modelica.Fluid.Types.Dynamics.SteadyState,
          massDynamics=if tau2 > Modelica.Constants.eps
                           then massDynamics else
                           Modelica.Fluid.Types.Dynamics.SteadyState,
          final p_start=p2_start,
          final T_start=T2_start,
          final X_start=X2_start,
          final C_start=C2_start,
          final C_nominal=C2_nominal) "Volume for fluid 2"
     annotation (Placement(transformation(
          origin={-66,-70},
          extent={{-10,-10},{10,10}},
          rotation=180)));

    Buildings.Fluid.FixedResistances.PressureDrop preDro1(
      redeclare final package Medium = Medium1,
      final m_flow_nominal=m1_flow_nominal,
      final deltaM=deltaM1,
      final allowFlowReversal=allowFlowReversal1,
      final show_T=false,
      final from_dp=from_dp1,
      final linearized=linearizeFlowResistance1,
      final homotopyInitialization=homotopyInitialization,
      final dp_nominal=dp1_nominal) "Flow resistance of fluid 1"
      annotation (Placement(transformation(extent={{-60,50},{-40,70}})));

    Buildings.Fluid.FixedResistances.PressureDrop preDro2(
      redeclare final package Medium = Medium2,
      final m_flow_nominal=m2_flow_nominal,
      final deltaM=deltaM2,
      final allowFlowReversal=allowFlowReversal2,
      final show_T=false,
      final from_dp=from_dp2,
      final linearized=linearizeFlowResistance2,
      final homotopyInitialization=homotopyInitialization,
      final dp_nominal=dp2_nominal) "Flow resistance of fluid 2"
      annotation (Placement(transformation(extent={{-26,-70},{-46,-50}})));

    BaseClasses.OutletConditions outletConditions
      annotation (Placement(transformation(extent={{0,20},{-20,0}})));
  public
    BaseClasses.Pad_Efficiency pad_Efficiency(tk_pad=PadThickness)
      annotation (Placement(transformation(extent={{42,-2},{22,18}})));
    Buildings.Fluid.Sensors.TemperatureWetBulbTwoPort senWetBul(
      redeclare package Medium = Medium2,
      allowFlowReversal=true,
      m_flow_nominal=m2_flow_nominal,
      tau=0,
      initType=Modelica.Blocks.Types.Init.SteadyState,
      TWetBul_start(displayUnit="K"))          annotation (Placement(transformation(extent={{44,-68},
              {28,-52}})));
    Buildings.Fluid.Sensors.TemperatureTwoPort senTem1(
      redeclare package Medium = Medium2,
      m_flow_nominal=m2_flow_nominal,
      tau=0,
      initType=Modelica.Blocks.Types.Init.InitialState,
      T_start(displayUnit="K"))                annotation (Placement(transformation(extent={{18,-66},
              {2,-50}})));
    Modelica.Blocks.Math.Product VolumeFlowRate2MassFlowRate annotation (Placement(transformation(extent={{4,-4},{
              -4,4}},
          rotation=90,
          origin={-52,-20})));
    parameter String substanceName "Name of species substance";
    Buildings.Fluid.Sensors.MassFractionTwoPort senMasFra(
      redeclare package Medium = Buildings.Media.Air,
      m_flow_nominal=m2_flow_nominal,
      tau=tau2,
      initType=Modelica.Blocks.Types.Init.NoInit)
                annotation (Placement(transformation(extent={{66,-68},{50,-52}})));
    BaseClasses.Water_consumption water_consumption(f_drift=DriftFactor, Rcon=
          Rcon)
      annotation (Placement(transformation(extent={{-60,20},{-80,0}})));
  protected
    parameter Medium1.ThermodynamicState sta1_nominal=Medium1.setState_pTX(
        T=Medium1.T_default, p=Medium1.p_default, X=Medium1.X_default);
    parameter Modelica.SIunits.Density rho1_nominal=Medium1.density(sta1_nominal)
      "Density, used to compute fluid volume";
    parameter Medium2.ThermodynamicState sta2_nominal=Medium2.setState_pTX(
        T=Medium2.T_default, p=Medium2.p_default, X=Medium2.X_default);
    parameter Modelica.SIunits.Density rho2_nominal=Medium2.density(sta2_nominal)
      "Density, used to compute fluid volume";

    parameter Medium1.ThermodynamicState sta1_start=Medium1.setState_pTX(
        T=T1_start, p=p1_start, X=X1_start);
    parameter Modelica.SIunits.SpecificEnthalpy h1_outflow_start = Medium1.specificEnthalpy(sta1_start)
      "Start value for outflowing enthalpy";
    parameter Medium2.ThermodynamicState sta2_start=Medium2.setState_pTX(
        T=T2_start, p=p2_start, X=X2_start);
    parameter Modelica.SIunits.SpecificEnthalpy h2_outflow_start = Medium2.specificEnthalpy(sta2_start)
      "Start value for outflowing enthalpy";

   Buildings.Fluid.Sensors.DensityTwoPort senDen(
      redeclare package Medium = Medium1,
      m_flow_nominal=m1_flow_nominal,
      tau=tau1,
      initType=Modelica.Blocks.Types.Init.InitialState)
     annotation (Placement(transformation(extent={{-30,68},{-14,52}})));
   Buildings.Fluid.Sensors.Velocity senVel(
      redeclare package Medium = Medium2,
      m_flow_nominal=m2_flow_nominal,
      tau=tau2,
      initType=Modelica.Blocks.Types.Init.InitialState,
      T_start(displayUnit="K"),
      A=PadArea)
     annotation (Placement(transformation(extent={{88,-68},{72,-52}})));
   Buildings.Fluid.Sensors.MassFlowRate senMasFlo(redeclare package Medium =
          Medium2)
     annotation (Placement(transformation(extent={{-4,-68},{-20,-52}})));
  initial equation
    // Check for tau1
    assert((energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) or
            tau1 > Modelica.Constants.eps,
  "The parameter tau1, or the volume of the model from which tau may be derived, is unreasonably small.
 You need to set energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState to model steady-state.
 Received tau1 = "   + String(tau1) + "\n");
    assert((massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) or
            tau1 > Modelica.Constants.eps,
  "The parameter tau1, or the volume of the model from which tau may be derived, is unreasonably small.
 You need to set massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState to model steady-state.
 Received tau1 = "   + String(tau1) + "\n");

   // Check for tau2
    assert((energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) or
            tau2 > Modelica.Constants.eps,
  "The parameter tau2, or the volume of the model from which tau may be derived, is unreasonably small.
 You need to set energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState to model steady-state.
 Received tau2 = "   + String(tau2) + "\n");
    assert((massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) or
            tau2 > Modelica.Constants.eps,
  "The parameter tau2, or the volume of the model from which tau may be derived, is unreasonably small.
 You need to set massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState to model steady-state.
 Received tau2 = "   + String(tau2) + "\n");

  equation
    connect(vol2.ports[1], port_b2) annotation (Line(
        points={{-64,-60},{-100,-60}},
        color={0,127,255}));
    connect(senDen.port_a, preDro1.port_b)
      annotation (Line(points={{-30,60},{-40,60}}, color={0,127,255}));
    connect(preDro2.port_b, vol2.ports[2])
      annotation (Line(points={{-46,-60},{-68,-60}}, color={0,127,255}));
    connect(preDro2.port_a, senMasFlo.port_b)
      annotation (Line(points={{-26,-60},{-20,-60}}, color={0,127,255}));
    connect(senMasFlo.port_a, senTem1.port_b)
      annotation (Line(points={{-4,-60},{-2,-60},{-2,-58},{2,-58}},
                                                  color={0,127,255}));
    connect(senTem1.port_a, senWetBul.port_b)
      annotation (Line(points={{18,-58},{24,-58},{24,-60},{28,-60}},
                                                   color={0,127,255}));
    connect(senVel.port_a, port_a2)
      annotation (Line(points={{88,-60},{100,-60}}, color={0,127,255}));
    connect(senVel.v, pad_Efficiency.vel_air) annotation (Line(points={{80,
            -51.2},{80,5.2},{41.4,5.2}},
                                  color={0,0,127}));
    connect(senWetBul.T, outletConditions.Twb_air) annotation (Line(points={{36,-51.2},
            {36,-6},{12,-6},{12,8.4},{-0.4,8.4}}, color={0,0,127}));
    connect(pad_Efficiency.efficiency, outletConditions.Efficiency) annotation (
        Line(points={{21.8,7.8},{16,7.8},{16,13.4},{-0.4,13.4}}, color={0,0,127}));
    connect(senTem1.T, outletConditions.Tdb_air) annotation (Line(points={{10,
            -49.2},{10,-10},{4,-10},{4,3.4},{-0.4,3.4}},
                                                  color={0,0,127}));
    connect(vol2.mWat_flow, VolumeFlowRate2MassFlowRate.y) annotation (Line(
          points={{-54,-78},{-54,-42},{-52,-42},{-52,-24.4}}, color={0,0,127}));
    connect(senDen.d, VolumeFlowRate2MassFlowRate.u1) annotation (Line(points={{-22,
            51.2},{-22,-6},{-54.4,-6},{-54.4,-15.2}},            color={0,0,127}));
    connect(vol1.ports[1], port_b1)
      annotation (Line(points={{34,60},{100,60}}, color={0,127,255}));
    connect(vol1.ports[2], senDen.port_b)
      annotation (Line(points={{38,60},{-14,60}}, color={0,127,255}));
    connect(senWetBul.port_a, senMasFra.port_b)
      annotation (Line(points={{44,-60},{50,-60}}, color={0,127,255}));
    connect(senMasFra.port_a, senVel.port_b)
      annotation (Line(points={{66,-60},{72,-60}}, color={0,127,255}));
    connect(port_a1, preDro1.port_a)
      annotation (Line(points={{-100,60},{-60,60}}, color={0,127,255}));
    connect(senDen.d, water_consumption.Density) annotation (Line(points={{-22,
            51.2},{-22,24},{-60.5,24},{-60.5,17.1}},color={0,0,127}));
    connect(water_consumption.Vol_evap, VolumeFlowRate2MassFlowRate.u2)
      annotation (Line(points={{-80,9.2},{-80,-2},{-49.6,-2},{-49.6,-15.2}},
          color={0,0,127}));
    connect(outletConditions.Wx_out, water_consumption.W_out) annotation (Line(
          points={{-20.4,9.2},{-32,9.2},{-32,12.5},{-60.5,12.5}}, color={0,0,
            127}));
    connect(senMasFra.X, water_consumption.W_in) annotation (Line(points={{58,
            -51.2},{58,-22},{-28,-22},{-28,3.3},{-60.5,3.3}}, color={0,0,127}));
    connect(senMasFlo.m_flow, water_consumption.m_flo) annotation (Line(points={{-12,
            -51.2},{-12,-30},{-38,-30},{-38,7.9},{-60.5,7.9}},       color={0,0,
            127}));
    annotation (
      Documentation(info="<html>
<p>
This component transports two fluid streams between four ports.
It provides the basic model for implementing a dynamic heat exchanger.
</p>
<p>
The model can be used as-is, although there will be no heat or mass transfer
between the two fluid streams.
To add heat transfer, heat flow can be added to the heat port of the two volumes.
See for example
<a href=\"Buildings.Fluid.Chillers.Carnot_y\">
Buildings.Fluid.Chillers.Carnot_y</a>.
To add moisture input into (or moisture output from) volume <code>vol2</code>,
the model can be replaced with
<a href=\"modelica://Buildings.Fluid.MixingVolumes.MixingVolumeMoistAir\">
Buildings.Fluid.MixingVolumes.MixingVolumeMoistAir</a>.
</p>
<h4>Implementation</h4>
<p>
The variable names follow the conventions used in
<a href=\"modelica://Modelica.Fluid.Examples.HeatExchanger.BaseClasses.BasicHX\">
Modelica.Fluid.Examples.HeatExchanger.BaseClasses.BasicHX</a>.
</p>
</html>",   revisions="<html>
<ul>
<li>
October 23, 2017, by Michael Wetter:<br/>
Made volume <code>vol1</code> replaceable. This is required for
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/1013\">Buildings, issue 1013</a>.
</li>
<li>
December 1, 2016, by Michael Wetter:<br/>
Updated model as <code>use_dh</code> is no longer a parameter in the pressure drop model.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/480\">#480</a>.
</li>
<li>
April 11, 2016 by Michael Wetter:<br/>
Corrected wrong hyperlink in documentation for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/450\">issue 450</a>.
</li>
<li>
January 26, 2016, by Michael Wetter:<br/>
Set <code>quantity</code> attributes.
</li>
<li>
November 13, 2015, by Michael Wetter:<br/>
Changed assignments of start values in <code>extends</code> statement.
This is for issue
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/299\">#299</a>.
</li>
<li>
June 2, 2015, by Filip Jorissen:<br/>
Removed final modifier from <code>mSenFac</code> in
<code>vol1</code> and <code>vol2</code>.
This is for issue
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/258=\">#258</a>.
</li>
<li>
May 6, 2015, by Michael Wetter:<br/>
Added missing propagation of <code>allowFlowReversal</code> to
instances <code>vol1</code> and <code>vol2</code>.
This is for issue
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/412\">#412</a>.
</li>
<li>
October 6, 2014, by Michael Wetter:<br/>
Changed medium declaration in pressure drop elements to be final.
</li>
<li>
May 28, 2014, by Michael Wetter:<br/>
Removed <code>annotation(Evaluate=true)</code> for parameters <code>tau1</code>
and <code>tau2</code>.
This is needed to allow changing the time constant after translation.
</li>
<li>
November 12, 2013, by Michael Wetter:<br/>
Removed <code>import Modelica.Constants</code> statement.
</li>
<li>
October 8, 2013, by Michael Wetter:<br/>
Removed parameter <code>show_V_flow</code>.
</li>
<li>
September 26, 2013, by Michael Wetter:<br/>
Removed unrequired <code>sum</code> operator.
</li>
<li>
February 6, 2012, by Michael Wetter:<br/>
Updated documentation.
</li>
<li>
February 3, 2012, by Michael Wetter:<br/>
Removed assignment of <code>m_flow_small</code> as it is no
longer used in its base class.
</li>
<li>
July 29, 2011, by Michael Wetter:
<ul>
<li>
Changed values of
<code>h_outflow_a1_start</code>,
<code>h_outflow_b1_start</code>,
<code>h_outflow_a2_start</code> and
<code>h_outflow_b2_start</code>, and
declared them as final.
</li>
<li>
Set nominal values for <code>vol1.C</code> and <code>vol2.C</code>.
</li>
</ul>
</li>
<li>
July 11, 2011, by Michael Wetter:<br/>
Changed parameterization of fluid volume so that steady-state balance is
used when <code>tau = 0</code>.
</li>
<li>
March 25, 2011, by Michael Wetter:<br/>
Added homotopy operator.
</li>
<li>
April 13, 2009, by Michael Wetter:<br/>
Added model to compute flow friction.
</li>
<li>
September 10, 2008 by Michael Wetter:<br/>
Added <code>stateSelect=StateSelect.always</code> for temperature of volume 1.
</li>
<li>
Changed temperature sensor from Celsius to Kelvin.
Unit conversion should be made during output
processing.
</li>
<li>
August 5, 2008, by Michael Wetter:<br/>
Replaced instances of <code>Delays.DelayFirstOrder</code> with instances of
<code>MixingVolumes.MixingVolume</code>. This allows to extract liquid for a condensing cooling
coil model.
</li>
<li>
March 25, 2008, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"),
      Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Rectangle(
            extent={{-70,80},{70,-80}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-99,64},{102,54}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-99,-56},{102,-66}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}));
  end DEC_EP_Method1;

  model DEC_Neumerical_Method2_test2

        extends Buildings.Fluid.Interfaces.PartialFourPortInterface(
        replaceable package Medium2 =
          Modelica.Media.Interfaces.PartialCondensingGases,
      allowFlowReversal2=false,
      allowFlowReversal1=false,
      port_a1(h_outflow(start=h1_outflow_start)),
      port_b1(h_outflow(start=h1_outflow_start)),
      port_a2(h_outflow(start=h2_outflow_start)),
      port_b2(h_outflow(start=h2_outflow_start)));
    extends Buildings.Fluid.Interfaces.FourPortFlowResistanceParameters(
       final computeFlowResistance1=true, final computeFlowResistance2=true);


  parameter Modelica.SIunits.Thickness Thickness=0.100
      "Evaporative cooling pad thickness(m)"
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

  parameter Modelica.SIunits.Thickness Length = 2.1   "Evaporative cooling pad length(m)"
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

  parameter Modelica.SIunits.Thickness Height = 0.91
                                                    "Evaporative cooling pad Height(m)"
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

   parameter Modelica.SIunits.Conductivity K_value= 0.04 " Evaporative cooling pad thermal conductivity"
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

   parameter Real Contact_surface_area = 440
                                            "Surface area per unit volume in comact with air(m2/m3)"
   annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

   parameter Real Pore_to_Volume_ratio = 0.93
                                             "Volume of pores to volume of the cooler pad (m3/m3)"
   annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

  parameter Modelica.SIunits.Area PadArea= Length * Height
      "Evaporative cooling pad area(m2)"
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

   parameter Real P_atm = 101325
                                "in pascals";

   parameter Real DriftFactor = 0.1
      "Drift factor[user input]"
    annotation (Dialog(group="Water consumption",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

   parameter Real Rcon = 2
      "Ratio of solids in the blowdown water[user input]"
    annotation (Dialog(group="Water consumption",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

    parameter Modelica.SIunits.Time tau1 = 1 "Time constant at nominal flow"
       annotation (Dialog(tab = "Dynamics", group="Nominal condition"));
    parameter Modelica.SIunits.Time tau2 = 1 "Time constant at nominal flow"
       annotation (Dialog(tab = "Dynamics", group="Nominal condition"));

    // Advanced
    parameter Boolean homotopyInitialization = true "= true, use homotopy method"
      annotation(Evaluate=true, Dialog(tab="Advanced"));

    // Assumptions
    parameter Modelica.Fluid.Types.Dynamics energyDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial
      "Type of energy balance: dynamic (3 initialization options) or steady state"
      annotation(Evaluate=true, Dialog(tab = "Dynamics", group="Equations"));
    parameter Modelica.Fluid.Types.Dynamics massDynamics=if tau1 > Modelica.Constants.eps
         then energyDynamics else Modelica.Fluid.Types.Dynamics.SteadyState
      "Type of mass balance: dynamic (3 initialization options) or steady state"
      annotation(Evaluate=true, Dialog(tab = "Dynamics", group="Equations"));

    // Initialization
    parameter Medium1.AbsolutePressure p1_start = Medium1.p_default
      "Start value of pressure"
      annotation(Dialog(tab = "Initialization", group = "Medium 1"));
    parameter Medium1.Temperature T1_start = Medium1.T_default
      "Start value of temperature"
      annotation(Dialog(tab = "Initialization", group = "Medium 1"));
    parameter Medium1.MassFraction X1_start[Medium1.nX] = Medium1.X_default
      "Start value of mass fractions m_i/m"
      annotation (Dialog(tab="Initialization", group = "Medium 1", enable=Medium1.nXi > 0));
    parameter Medium1.ExtraProperty C1_start[Medium1.nC](
      final quantity=Medium1.extraPropertiesNames)=fill(0, Medium1.nC)
      "Start value of trace substances"
      annotation (Dialog(tab="Initialization", group = "Medium 1", enable=Medium1.nC > 0));
    parameter Medium1.ExtraProperty C1_nominal[Medium1.nC](
      final quantity=Medium1.extraPropertiesNames) = fill(1E-2, Medium1.nC)
      "Nominal value of trace substances. (Set to typical order of magnitude.)"
     annotation (Dialog(tab="Initialization", group = "Medium 1", enable=Medium1.nC > 0));

    parameter Medium2.AbsolutePressure p2_start = Medium2.p_default
      "Start value of pressure"
      annotation(Dialog(tab = "Initialization", group = "Medium 2"));
    parameter Medium2.Temperature T2_start = Medium2.T_default
      "Start value of temperature"
      annotation(Dialog(tab = "Initialization", group = "Medium 2"));
    parameter Medium2.MassFraction X2_start[Medium2.nX] = Medium2.X_default
      "Start value of mass fractions m_i/m"
      annotation (Dialog(tab="Initialization", group = "Medium 2", enable=Medium2.nXi > 0));
    parameter Medium2.ExtraProperty C2_start[Medium2.nC](
      final quantity=Medium2.extraPropertiesNames)=fill(0, Medium2.nC)
      "Start value of trace substances"
      annotation (Dialog(tab="Initialization", group = "Medium 2", enable=Medium2.nC > 0));
    parameter Medium2.ExtraProperty C2_nominal[Medium2.nC](
      final quantity=Medium2.extraPropertiesNames) = fill(1E-2, Medium2.nC)
      "Nominal value of trace substances. (Set to typical order of magnitude.)"
     annotation (Dialog(tab="Initialization", group = "Medium 2", enable=Medium2.nC > 0));

    Modelica.SIunits.HeatFlowRate Q1_flow=volume_1_water.heatPort.Q_flow
      "Heat flow rate into medium 1";
    Modelica.SIunits.HeatFlowRate Q2_flow=Volume_2_air.heatPort.Q_flow
      "Heat flow rate into medium 2";

    replaceable
      Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatMoisturePort volume_1_water(nPorts=2)
      constrainedby
      Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatPort(
      redeclare final package Medium = Medium1,
      nPorts=2,
      V=m1_flow_nominal*tau1/rho1_nominal,
      final allowFlowReversal=allowFlowReversal1,
      final m_flow_nominal=m1_flow_nominal,
      energyDynamics=if tau1 > Modelica.Constants.eps then energyDynamics else
          Modelica.Fluid.Types.Dynamics.SteadyState,
      massDynamics=if tau1 > Modelica.Constants.eps then massDynamics else
          Modelica.Fluid.Types.Dynamics.SteadyState,
      final p_start=p1_start,
      final T_start=T1_start,
      final X_start=X1_start,
      final C_start=C1_start,
      final C_nominal=C1_nominal,
      mSenFac=1) "Volume for fluid 1_water"
      annotation (Placement(transformation(extent={{12,60},{30,42}})));

    replaceable
      Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatMoisturePort Volume_2_air(nPorts=2)
      constrainedby
      Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatPort(
      redeclare final package Medium = Medium2,
      nPorts=2,
      V=m2_flow_nominal*tau2/rho2_nominal,
      final allowFlowReversal=allowFlowReversal2,
      mSenFac=1,
      final m_flow_nominal=m2_flow_nominal,
      energyDynamics=if tau2 > Modelica.Constants.eps then energyDynamics else
          Modelica.Fluid.Types.Dynamics.SteadyState,
      massDynamics=if tau2 > Modelica.Constants.eps then massDynamics else
          Modelica.Fluid.Types.Dynamics.SteadyState,
      final p_start=p2_start,
      final T_start=T2_start,
      final X_start=X2_start,
      final C_start=C2_start,
      final C_nominal=C2_nominal) "Volume for fluid 2_air" annotation (Placement(
          transformation(
          origin={-50,-70},
          extent={{-10,-10},{10,10}},
          rotation=180)));

    Buildings.Fluid.FixedResistances.PressureDrop preDro1(
      redeclare final package Medium = Medium1,
      final m_flow_nominal=m1_flow_nominal,
      final deltaM=deltaM1,
      final allowFlowReversal=allowFlowReversal1,
      final show_T=false,
      final from_dp=from_dp1,
      final linearized=linearizeFlowResistance1,
      final homotopyInitialization=homotopyInitialization,
      final dp_nominal=dp1_nominal) "Flow resistance of water side"
      annotation (Placement(transformation(extent={{42,54},{54,66}})));

    Buildings.Fluid.FixedResistances.PressureDrop preDro2(
      redeclare final package Medium = Medium2,
      final m_flow_nominal=m2_flow_nominal,
      final deltaM=deltaM2,
      final allowFlowReversal=allowFlowReversal2,
      final show_T=false,
      final from_dp=from_dp2,
      final linearized=linearizeFlowResistance2,
      final homotopyInitialization=homotopyInitialization,
      final dp_nominal=dp2_nominal) "Flow resistance of air side"
      annotation (Placement(transformation(extent={{-74,-66},{-86,-54}})));

    Buildings.Fluid.Sensors.TemperatureWetBulbTwoPort senWetBul2(
      redeclare package Medium = Medium2,
      allowFlowReversal=true,
      m_flow_nominal=m2_flow_nominal,
      tau=0,
      initType=Modelica.Blocks.Types.Init.SteadyState,
      TWetBul_start(displayUnit="K")) "Air side WBT (~ambient WBT)"
      annotation (Placement(transformation(extent={{76,-68},{60,-52}})));
    Buildings.Fluid.Sensors.TemperatureTwoPort senTem2(
      redeclare package Medium = Medium2,
      m_flow_nominal=m2_flow_nominal,
      tau=0,
      initType=Modelica.Blocks.Types.Init.InitialState,
      T_start(displayUnit="K")) "Air side temperature (~ ambient DBT)"
      annotation (Placement(transformation(extent={{94,-68},{78,-52}})));
    Buildings.Fluid.Sensors.MassFractionTwoPort senMasFra2(
      redeclare package Medium = Buildings.Media.Air,
      m_flow_nominal=m2_flow_nominal,
      tau=0,
      initType=Modelica.Blocks.Types.Init.NoInit)
      "Air side inlet humidity ratio kg/kg"
      annotation (Placement(transformation(extent={{4,-66},{-8,-54}})));
    BaseClasses.HeatTransfer heatTransfer(Area_s=PadArea)
      annotation (Placement(transformation(extent={{12,-40},{-10,-18}})));
    Buildings.Fluid.Sensors.TemperatureTwoPort senTem1(
      redeclare package Medium = Medium1,
      m_flow_nominal=m1_flow_nominal,
      tau=0,
      initType=Modelica.Blocks.Types.Init.InitialState,
      T_start(displayUnit="K")) "Water side fluid temperature"
                                               annotation (Placement(transformation(extent={{74,68},
              {60,52}})));
    BaseClasses.NTU_Eff_calculation nTU_Eff_calculation(
      p_atm=P_atm,
      k=K_value,
      Thickness=Thickness,
      length=Length,
      height=Height,
      p_V=Contact_surface_area)
      annotation (Placement(transformation(extent={{56,14},{30,-12}})));
  protected
    parameter Medium1.ThermodynamicState sta1_nominal=Medium1.setState_pTX(
        T=Medium1.T_default, p=Medium1.p_default, X=Medium1.X_default);
    parameter Modelica.SIunits.Density rho1_nominal=Medium1.density(sta1_nominal)
      "Density, used to compute fluid volume";
    parameter Medium2.ThermodynamicState sta2_nominal=Medium2.setState_pTX(
        T=Medium2.T_default, p=Medium2.p_default, X=Medium2.X_default);
    parameter Modelica.SIunits.Density rho2_nominal=Medium2.density(sta2_nominal)
      "Density, used to compute fluid volume";

    parameter Medium1.ThermodynamicState sta1_start=Medium1.setState_pTX(
        T=T1_start, p=p1_start, X=X1_start);
    parameter Modelica.SIunits.SpecificEnthalpy h1_outflow_start = Medium1.specificEnthalpy(sta1_start)
      "Start value for outflowing enthalpy";
    parameter Medium2.ThermodynamicState sta2_start=Medium2.setState_pTX(
        T=T2_start, p=p2_start, X=X2_start);
    parameter Modelica.SIunits.SpecificEnthalpy h2_outflow_start = Medium2.specificEnthalpy(sta2_start)
      "Start value for outflowing enthalpy";

   Buildings.Fluid.Sensors.DensityTwoPort senDen1(
      redeclare package Medium = Medium1,
      m_flow_nominal=m1_flow_nominal,
      tau=tau1,
      initType=Modelica.Blocks.Types.Init.InitialState) "Water side density"
      annotation (Placement(transformation(extent={{-58,66},{-46,54}})));
   Buildings.Fluid.Sensors.Velocity senVel2(
      redeclare package Medium = Medium2,
      m_flow_nominal=m2_flow_nominal,
      tau=0,
      initType=Modelica.Blocks.Types.Init.InitialState,
      A=PadArea) "Air side velocity over the cooling pad"
      annotation (Placement(transformation(extent={{40,-66},{28,-54}})));
   Buildings.Fluid.Sensors.MassFlowRate senMasFlo2(redeclare package Medium =
          Medium2) "Air side mass flow rate"
      annotation (Placement(transformation(extent={{56,-66},{44,-54}})));
   Buildings.Fluid.Sensors.MassFlowRate senMasFlo1(redeclare package Medium =
          Buildings.Media.Water) "Water side mass flow rate"
     annotation (Placement(transformation(extent={{6,-6},{-6,6}},
          rotation=180,
          origin={-6,60})));
  public
    BaseClasses.Water_consumption water_consumption(f_drift=DriftFactor, Rcon=
          Rcon)
      annotation (Placement(transformation(extent={{-20,0},{-40,-20}})));
    Modelica.Blocks.Math.Product VolumeFlowRate2MassFlowRate annotation (Placement(transformation(extent={{-4,-4},
              {4,4}},
          rotation=180,
          origin={-62,-16})));
    Modelica.Blocks.Math.Product VolumeFlowRate2MassFlowRate1
                                                             annotation (Placement(transformation(extent={{4,-4},{
              -4,4}},
          rotation=0,
          origin={-62,-2})));
    Modelica.Blocks.Math.Gain gain(k=-1)
      annotation (Placement(transformation(extent={{-74,-4},{-80,2}})));
  initial equation
    // Check for tau1
    assert((energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) or
            tau1 > Modelica.Constants.eps,
  "The parameter tau1, or the volume of the model from which tau may be derived, is unreasonably small.
 You need to set energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState to model steady-state.
 Received tau1 = "   + String(tau1) + "\n");
    assert((massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) or
            tau1 > Modelica.Constants.eps,
  "The parameter tau1, or the volume of the model from which tau may be derived, is unreasonably small.
 You need to set massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState to model steady-state.
 Received tau1 = "   + String(tau1) + "\n");

   // Check for tau2
    assert((energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) or
            tau2 > Modelica.Constants.eps,
  "The parameter tau2, or the volume of the model from which tau may be derived, is unreasonably small.
 You need to set energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState to model steady-state.
 Received tau2 = "   + String(tau2) + "\n");
    assert((massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) or
            tau2 > Modelica.Constants.eps,
  "The parameter tau2, or the volume of the model from which tau may be derived, is unreasonably small.
 You need to set massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState to model steady-state.
 Received tau2 = "   + String(tau2) + "\n");

  equation

    connect(senDen1.port_b, senMasFlo1.port_a)
      annotation (Line(points={{-46,60},{-12,60}}, color={0,127,255}));
    connect(senMasFlo1.port_b, volume_1_water.ports[1])
      annotation (Line(points={{-8.88178e-16,60},{19.2,60}}, color={0,127,255}));
    connect(gain.y, volume_1_water.mWat_flow) annotation (Line(points={{-80.3,-1},
            {-84,-1},{-84,44},{-36,44},{-36,43.8},{10.2,43.8}}, color={0,0,127}));
    connect(port_b1, senTem1.port_a)
      annotation (Line(points={{100,60},{74,60}}, color={0,127,255}));
    connect(senTem1.T, nTU_Eff_calculation.Tw1) annotation (Line(points={{67,51.2},
            {67,7.37},{55.61,7.37}}, color={0,0,127}));
    connect(senTem1.port_b, preDro1.port_b)
      annotation (Line(points={{60,60},{54,60}}, color={0,127,255}));
    connect(volume_1_water.ports[2], preDro1.port_a)
      annotation (Line(points={{22.8,60},{42,60}}, color={0,127,255}));
    connect(senVel2.v, nTU_Eff_calculation.v_a) annotation (Line(points={{34,-53.4},
            {34,-34.7},{47.55,-34.7},{47.55,-12.13}}, color={0,0,127}));
    connect(port_a2, senTem2.port_a)
      annotation (Line(points={{100,-60},{94,-60}}, color={0,127,255}));
    connect(senTem2.port_b, senWetBul2.port_a)
      annotation (Line(points={{78,-60},{76,-60}}, color={0,127,255}));
    connect(senWetBul2.T, nTU_Eff_calculation.Ta_wb1) annotation (Line(points={{68,
            -51.2},{68,-1.73},{55.61,-1.73}}, color={0,0,127}));
    connect(senTem2.T, nTU_Eff_calculation.Ta_db1) annotation (Line(points={{86,-51.2},
            {82,-51.2},{82,2.69},{55.61,2.69}}, color={0,0,127}));
    connect(senWetBul2.port_b, senMasFlo2.port_a)
      annotation (Line(points={{60,-60},{56,-60}}, color={0,127,255}));
    connect(senVel2.port_a, senMasFlo2.port_b)
      annotation (Line(points={{40,-60},{44,-60}}, color={0,127,255}));
    connect(senMasFlo2.m_flow, nTU_Eff_calculation.ma_flow) annotation (Line(
          points={{50,-53.4},{50,-42},{60,-42},{60,-6.15},{55.61,-6.15}}, color={0,
            0,127}));
    connect(senVel2.port_b, senMasFra2.port_a)
      annotation (Line(points={{28,-60},{4,-60}}, color={0,127,255}));
    connect(senMasFra2.X, water_consumption.W_in) annotation (Line(points={{-2,-53.4},
            {-14,-53.4},{-14,-16.7},{-20.5,-16.7}}, color={0,0,127}));
    connect(nTU_Eff_calculation.X_2, water_consumption.W_out) annotation (Line(
          points={{29.74,-7.84},{-20.5,-7.84},{-20.5,-7.5}}, color={0,0,127}));
    connect(preDro2.port_b, port_b2)
      annotation (Line(points={{-86,-60},{-100,-60}}, color={0,127,255}));
    connect(senDen1.d, VolumeFlowRate2MassFlowRate1.u1) annotation (Line(points={{
            -52,53.4},{-54,53.4},{-54,0.4},{-57.2,0.4}}, color={0,0,127}));
    connect(senDen1.d, water_consumption.Density) annotation (Line(points={{-52,53.4},
            {-54,53.4},{-54,4},{-14,4},{-14,-2.9},{-20.5,-2.9}}, color={0,0,127}));
    connect(VolumeFlowRate2MassFlowRate1.y, gain.u) annotation (Line(points={{-66.4,
            -2},{-70,-2},{-70,-1},{-73.4,-1}}, color={0,0,127}));
    connect(water_consumption.Vol_water, VolumeFlowRate2MassFlowRate1.u2)
      annotation (Line(points={{-40.2,-6.8},{-44,-6.8},{-44,-4.4},{-57.2,-4.4}},
          color={0,0,127}));
    connect(water_consumption.Vol_evap, VolumeFlowRate2MassFlowRate.u1)
      annotation (Line(points={{-40,-10.8},{-44,-10.8},{-44,-18.4},{-57.2,-18.4}},
          color={0,0,127}));
    connect(VolumeFlowRate2MassFlowRate.y, Volume_2_air.mWat_flow) annotation (
        Line(points={{-66.4,-16},{-68,-16},{-68,-78},{-38,-78}}, color={0,0,127}));
    connect(port_a1, senDen1.port_a)
      annotation (Line(points={{-100,60},{-58,60}}, color={0,127,255}));
    connect(senMasFlo1.m_flow, nTU_Eff_calculation.mw_flow) annotation (Line(
          points={{-6,53.4},{-8,53.4},{-8,20},{47.03,20},{47.03,13.35}}, color={0,
            0,127}));
    connect(nTU_Eff_calculation.hc, heatTransfer.hc) annotation (Line(points={{29.74,
            8.28},{22,8.28},{22,-18.4583},{13.7,-18.4583}}, color={0,0,127}));
    connect(nTU_Eff_calculation.hm, heatTransfer.hm) annotation (Line(points={{29.74,
            -1.08},{24,-1.08},{24,-30.8333},{13.6,-30.8333}}, color={0,0,127}));
    connect(nTU_Eff_calculation.Ta_db2, heatTransfer.T_db2) annotation (Line(
          points={{29.74,-4.46},{26,-4.46},{26,-34.1333},{13.6,-34.1333}}, color={
            0,0,127}));
    connect(nTU_Eff_calculation.X_2, heatTransfer.Wx_2) annotation (Line(points={{29.74,
            -7.84},{29.74,-37.6167},{13.6,-37.6167}},   color={0,0,127}));
    connect(senMasFra2.X, heatTransfer.Wx_1) annotation (Line(points={{-2,-53.4},
            {10,-53.4},{10,-52},{22,-52},{22,-24.9667},{13.6,-24.9667}},
                                                                  color={0,0,127}));
    connect(senTem2.T, heatTransfer.T_db1) annotation (Line(points={{86,-51.2},
            {86,-52},{82,-52},{82,-21.6667},{13.6,-21.6667}},
                                                    color={0,0,127}));
    connect(preDro2.port_a, Volume_2_air.ports[1])
      annotation (Line(points={{-74,-60},{-48,-60}}, color={0,127,255}));
    connect(Volume_2_air.ports[2], senMasFra2.port_b)
      annotation (Line(points={{-52,-60},{-8,-60}}, color={0,127,255}));
    connect(senDen1.d, VolumeFlowRate2MassFlowRate.u2) annotation (Line(points=
            {{-52,53.4},{-52,-13.6},{-57.2,-13.6}}, color={0,0,127}));
    connect(senMasFlo2.m_flow, water_consumption.m_flo) annotation (Line(points=
           {{50,-53.4},{48,-53.4},{48,-42},{-12,-42},{-12,-12.1},{-20.5,-12.1}},
          color={0,0,127}));
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)),
                Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)),
               Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)),
              Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)),
              Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)),
                              Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)),
      Documentation(info="<html>
<p>
This component transports two fluid streams between four ports.
It provides the basic model for implementing a dynamic heat exchanger.
</p>
<p>
The model can be used as-is, although there will be no heat or mass transfer
between the two fluid streams.
To add heat transfer, heat flow can be added to the heat port of the two volumes.
See for example
<a href=\"Buildings.Fluid.Chillers.Carnot_y\">
Buildings.Fluid.Chillers.Carnot_y</a>.
To add moisture input into (or moisture output from) volume <code>vol2</code>,
the model can be replaced with
<a href=\"modelica://Buildings.Fluid.MixingVolumes.MixingVolumeMoistAir\">
Buildings.Fluid.MixingVolumes.MixingVolumeMoistAir</a>.
</p>
<h4>Implementation</h4>
<p>
The variable names follow the conventions used in
<a href=\"modelica://Modelica.Fluid.Examples.HeatExchanger.BaseClasses.BasicHX\">
Modelica.Fluid.Examples.HeatExchanger.BaseClasses.BasicHX</a>.
</p>
</html>",   revisions="<html>
<ul>
<li>
October 23, 2017, by Michael Wetter:<br/>
Made volume <code>vol1</code> replaceable. This is required for
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/1013\">Buildings, issue 1013</a>.
</li>
<li>
December 1, 2016, by Michael Wetter:<br/>
Updated model as <code>use_dh</code> is no longer a parameter in the pressure drop model.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/480\">#480</a>.
</li>
<li>
April 11, 2016 by Michael Wetter:<br/>
Corrected wrong hyperlink in documentation for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/450\">issue 450</a>.
</li>
<li>
January 26, 2016, by Michael Wetter:<br/>
Set <code>quantity</code> attributes.
</li>
<li>
November 13, 2015, by Michael Wetter:<br/>
Changed assignments of start values in <code>extends</code> statement.
This is for issue
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/299\">#299</a>.
</li>
<li>
June 2, 2015, by Filip Jorissen:<br/>
Removed final modifier from <code>mSenFac</code> in
<code>vol1</code> and <code>vol2</code>.
This is for issue
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/258=\">#258</a>.
</li>
<li>
May 6, 2015, by Michael Wetter:<br/>
Added missing propagation of <code>allowFlowReversal</code> to
instances <code>vol1</code> and <code>vol2</code>.
This is for issue
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/412\">#412</a>.
</li>
<li>
October 6, 2014, by Michael Wetter:<br/>
Changed medium declaration in pressure drop elements to be final.
</li>
<li>
May 28, 2014, by Michael Wetter:<br/>
Removed <code>annotation(Evaluate=true)</code> for parameters <code>tau1</code>
and <code>tau2</code>.
This is needed to allow changing the time constant after translation.
</li>
<li>
November 12, 2013, by Michael Wetter:<br/>
Removed <code>import Modelica.Constants</code> statement.
</li>
<li>
October 8, 2013, by Michael Wetter:<br/>
Removed parameter <code>show_V_flow</code>.
</li>
<li>
September 26, 2013, by Michael Wetter:<br/>
Removed unrequired <code>sum</code> operator.
</li>
<li>
February 6, 2012, by Michael Wetter:<br/>
Updated documentation.
</li>
<li>
February 3, 2012, by Michael Wetter:<br/>
Removed assignment of <code>m_flow_small</code> as it is no
longer used in its base class.
</li>
<li>
July 29, 2011, by Michael Wetter:
<ul>
<li>
Changed values of
<code>h_outflow_a1_start</code>,
<code>h_outflow_b1_start</code>,
<code>h_outflow_a2_start</code> and
<code>h_outflow_b2_start</code>, and
declared them as final.
</li>
<li>
Set nominal values for <code>vol1.C</code> and <code>vol2.C</code>.
</li>
</ul>
</li>
<li>
July 11, 2011, by Michael Wetter:<br/>
Changed parameterization of fluid volume so that steady-state balance is
used when <code>tau = 0</code>.
</li>
<li>
March 25, 2011, by Michael Wetter:<br/>
Added homotopy operator.
</li>
<li>
April 13, 2009, by Michael Wetter:<br/>
Added model to compute flow friction.
</li>
<li>
September 10, 2008 by Michael Wetter:<br/>
Added <code>stateSelect=StateSelect.always</code> for temperature of volume 1.
</li>
<li>
Changed temperature sensor from Celsius to Kelvin.
Unit conversion should be made during output
processing.
</li>
<li>
August 5, 2008, by Michael Wetter:<br/>
Replaced instances of <code>Delays.DelayFirstOrder</code> with instances of
<code>MixingVolumes.MixingVolume</code>. This allows to extract liquid for a condensing cooling
coil model.
</li>
<li>
March 25, 2008, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"),
      Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Rectangle(
            extent={{-70,80},{70,-80}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-99,64},{102,54}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-99,-56},{102,-66}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}),
      Diagram(graphics={Text(
            extent={{-96,86},{-70,76}},
            lineColor={28,108,200},
            textString="Water Side"), Text(
            extent={{-90,-86},{-68,-94}},
            lineColor={28,108,200},
            textString="Air Side")}));
  end DEC_Neumerical_Method2_test2;

  model DirectEvaporativeCooler
   extends EnergyPlus_DEC_Blocks.BaseClasses.FourPortHeatMassExchanger(
     redeclare final Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatPort vol1(
       final energyDynamics=energyDynamics,
       final massDynamics=energyDynamics,
       prescribedHeatFlowRate=false),
     redeclare final Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatMoisturePort vol2(
       final energyDynamics=energyDynamics,
       final massDynamics=energyDynamics,
       prescribedHeatFlowRate=false),
     replaceable package Medium1 = Buildings.Media.Water,
     replaceable package Medium2 = Buildings.Media.Air);


  parameter Modelica.SIunits.Thickness Thickness=0.100
      "Evaporative cooling pad thickness(m)"
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

  parameter Modelica.SIunits.Thickness Length = 2.1   "Evaporative cooling pad length(m)"
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

  parameter Modelica.SIunits.Thickness Height = 0.91
                                                    "Evaporative cooling pad Height(m)"
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

   parameter Modelica.SIunits.Conductivity K_value= 0.04 " Evaporative cooling pad thermal conductivity"
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

   parameter Real Contact_surface_area = 440
                                            "Surface area per unit volume in comact with air(m2/m3)"
   annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

   parameter Real Pore_to_Volume_ratio = 0.93
                                             "Volume of pores to volume of the cooler pad (m3/m3)"
   annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

  parameter Modelica.SIunits.Area PadArea= Length * Height
      "Evaporative cooling pad area(m2)"
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

   parameter Real P_atm = 101325
                                "in pascals";

   parameter Real DriftFactor = 0.1
      "Drift factor[user input]"
    annotation (Dialog(group="Water consumption",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

   parameter Real Rcon = 2
      "Ratio of solids in the blowdown water[user input]";


    BaseClasses.NTU_Eff_calculation nTU_Eff_calculation(
      p_atm=P_atm,
      k=K_value,
      Thickness=Thickness,
      length=Length,
      height=Height,
      p_V=Contact_surface_area)
      annotation (Placement(transformation(extent={{60,-40},{40,-20}})));

  protected
    Modelica.Blocks.Math.Product pro
      "Product to compute the latent heat flow rate"
      annotation (Placement(transformation(extent={{70,4},{50,24}})));
    Modelica.Blocks.Sources.RealExpression h_fg(final y=PadArea*Buildings.Utilities.Psychrometrics.Constants.h_fg)
      "Enthalpy of vaporization"
      annotation (Placement(transformation(extent={{100,10},{80,30}})));
    Buildings.HeatTransfer.Sources.PrescribedHeatFlow heaConVapAir
      "Heat conductor for latent heat flow rate, accounting for latent heat removed with vapor"
      annotation (Placement(transformation(extent={{0,-40},{-20,-20}})));
    Buildings.HeatTransfer.Sources.PrescribedHeatFlow convHea
      "Convective heat transfer"
      annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
    Modelica.Blocks.Math.Product proConv
      "Product to compute the convective heat flow rate"
      annotation (Placement(transformation(extent={{0,0},{-20,20}})));
    Modelica.Blocks.Sources.RealExpression dTdA(final y=PadArea*(senTem_b2.T -
          senTem_a2.T)) "Temperature difference"
      annotation (Placement(transformation(extent={{18,-24},{-2,-4}})));

  equation

    connect(nTU_Eff_calculation.hm, vol2.mWat_flow) annotation (Line(points={{39.8,
            -28.4},{30,-28.4},{30,-52},{-28,-52}},color={0,0,127}));
    connect(pro.u1,h_fg. y)
      annotation (Line(points={{72,20},{79,20}}, color={0,0,127}));
    connect(nTU_Eff_calculation.hm, pro.u2) annotation (Line(points={{39.8,-28.4},
            {30,-28.4},{30,-6},{80,-6},{80,8},{72,8}}, color={0,0,127}));
    connect(heaConVapAir.Q_flow, pro.y) annotation (Line(points={{0,-30},{24,-30},
            {24,14},{49,14}}, color={0,0,127}));
    connect(heaConVapAir.port, vol2.heatPort) annotation (Line(points={{-20,-30},{
            -32,-30},{-32,-40},{20,-40},{20,-60},{-30,-60}},color={191,0,0}));
    connect(nTU_Eff_calculation.hc, proConv.u1) annotation (Line(points={{39.8,-35.6},
            {28,-35.6},{28,16},{2,16}}, color={0,0,127}));
    connect(dTdA.y, proConv.u2) annotation (Line(points={{-3,-14},{-16,-14},{-16,-2},
            {14,-2},{14,4},{2,4}}, color={0,0,127}));
    connect(proConv.y, convHea.Q_flow) annotation (Line(points={{-21,10},{-90,10},
            {-90,-20},{-80,-20}}, color={0,0,127}));
    connect(convHea.port, vol2.heatPort) annotation (Line(points={{-60,-20},{-40,-20},
            {-40,-44},{22,-44},{22,-60},{-30,-60}}, color={191,0,0}));
    connect(senMasFlo1.m_flow, nTU_Eff_calculation.mw_flow) annotation (Line(
          points={{-74,91},{-74,98},{-60,98},{-60,40},{46,40},{46,-48},{53.1,-48},
            {53.1,-39.5}}, color={0,0,127}));
    connect(senMasFlo2.m_flow, nTU_Eff_calculation.ma_flow) annotation (Line(
          points={{80,-69},{72,-69},{72,-24.5},{59.7,-24.5}}, color={0,0,127}));
    connect(senVel.v, nTU_Eff_calculation.v_a) annotation (Line(points={{54,-69},{
            54,-64},{70,-64},{70,-14},{53.5,-14},{53.5,-19.9}}, color={0,0,127}));
    connect(senWetBul.T, nTU_Eff_calculation.Ta_wb1) annotation (Line(points={{28,
            -69},{28,-60},{68,-60},{68,-27.9},{59.7,-27.9}}, color={0,0,127}));
    connect(senTem_a2.T, nTU_Eff_calculation.Ta_db1) annotation (Line(points={{87,
            -50},{96,-50},{96,-31.3},{59.7,-31.3}}, color={0,0,127}));
    connect(senTem_a1.T, nTU_Eff_calculation.Tw1) annotation (Line(points={{-81,90},
            {-76,90},{-76,100},{78,100},{78,-34.9},{59.7,-34.9}}, color={0,0,127}));
    annotation (Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)),
                Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)),
               Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)),
              Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)),
              Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)),
                              Dialog(group="Evaporative cooler pad",
                  enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)),
      Documentation(info="",
            revisions=""),
      Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Rectangle(
            extent={{-70,80},{70,-80}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-99,64},{102,54}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-99,-56},{102,-66}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}),
      Diagram(graphics={Text(
            extent={{-62,98},{-36,88}},
            lineColor={28,108,200},
            textString="Water Side"), Text(
            extent={{-90,-86},{-68,-94}},
            lineColor={28,108,200},
            textString="Air Side")}));
  end DirectEvaporativeCooler;
  annotation (uses(Modelica(version="3.2.2"), Buildings(version="5.0.0")),
    version="1",
    conversion(noneFromVersion=""));
end EnergyPlus_DEC_Blocks;
