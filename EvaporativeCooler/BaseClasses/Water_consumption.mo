within EvaporativeCooler.BaseClasses;
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
