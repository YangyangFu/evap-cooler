within EvaporativeCooler.BaseClasses;
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
  Nu=  0.192*((l/Thickness)^(-0.191))*((max(1e-06,(Ta_db1- Ta_wb1))/Tw1)^(-0.03))* (Re^0.682) * (Pr^0.33);
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
