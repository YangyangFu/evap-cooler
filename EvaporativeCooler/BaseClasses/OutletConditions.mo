within EvaporativeCooler.BaseClasses;
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
