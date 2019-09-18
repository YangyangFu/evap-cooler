within EvaporativeCooler.BaseClasses;
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
