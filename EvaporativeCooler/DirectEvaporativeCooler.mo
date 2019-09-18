within EvaporativeCooler;
model DirectEvaporativeCooler
 extends EvaporativeCooler.BaseClasses.FourPortHeatMassExchanger(
    redeclare final
      Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatPort vol1(
      final energyDynamics=energyDynamics,
      final massDynamics=energyDynamics,
      prescribedHeatFlowRate=false),
    redeclare final
      Buildings.Fluid.MixingVolumes.BaseClasses.MixingVolumeHeatMoisturePort
      vol2(
      final energyDynamics=energyDynamics,
      final massDynamics=energyDynamics,
      prescribedHeatFlowRate=false),
    replaceable package Medium1 = Buildings.Media.Water,
    replaceable package Medium2 = Buildings.Media.Air,
    final A = PadArea);

parameter Modelica.SIunits.Thickness Thickness=0.100
    "Evaporative cooling pad thickness(m)"
  annotation (Dialog(group="Evaporative cooler pad",
                enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

parameter Modelica.SIunits.Thickness Length = 2.1   "Evaporative cooling pad length(m)"
  annotation (Dialog(group="Evaporative cooler pad",
                enable=not (energyDynamics==Modelica.Fluid.Types.Dynamics.SteadyState)));

parameter Modelica.SIunits.Thickness Height = 0.91 "Evaporative cooling pad Height(m)"
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

  BaseClasses.NTU_Eff_calculation eff(
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
  Modelica.Blocks.Sources.RealExpression h_fg(final y=PadArea*Buildings.Utilities.Psychrometrics.Constants.h_fg
        *(eff.X_2 - senMasFra.X))
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

  connect(eff.hm, vol2.mWat_flow) annotation (Line(points={{39.8,-28.4},{30,-28.4},
          {30,-52},{-28,-52}}, color={0,0,127}));
  connect(pro.u1,h_fg. y)
    annotation (Line(points={{72,20},{79,20}}, color={0,0,127}));
  connect(eff.hm, pro.u2) annotation (Line(points={{39.8,-28.4},{30,-28.4},{30,
          -6},{80,-6},{80,8},{72,8}}, color={0,0,127}));
  connect(heaConVapAir.Q_flow, pro.y) annotation (Line(points={{0,-30},{24,-30},
          {24,14},{49,14}}, color={0,0,127}));
  connect(heaConVapAir.port, vol2.heatPort) annotation (Line(points={{-20,-30},{
          -32,-30},{-32,-40},{20,-40},{20,-60},{-30,-60}},color={191,0,0}));
  connect(eff.hc, proConv.u1) annotation (Line(points={{39.8,-35.6},{28,-35.6},
          {28,16},{2,16}}, color={0,0,127}));
  connect(dTdA.y, proConv.u2) annotation (Line(points={{-3,-14},{-16,-14},{-16,-2},
          {14,-2},{14,4},{2,4}}, color={0,0,127}));
  connect(proConv.y, convHea.Q_flow) annotation (Line(points={{-21,10},{-90,10},
          {-90,-20},{-80,-20}}, color={0,0,127}));
  connect(convHea.port, vol2.heatPort) annotation (Line(points={{-60,-20},{-40,-20},
          {-40,-44},{22,-44},{22,-60},{-30,-60}}, color={191,0,0}));
  connect(senMasFlo1.m_flow, eff.mw_flow) annotation (Line(points={{-74,91},{-74,
          98},{-60,98},{-60,40},{46,40},{46,-48},{53.1,-48},{53.1,-39.5}},
        color={0,0,127}));
  connect(senMasFlo2.m_flow, eff.ma_flow) annotation (Line(points={{80,-69},{72,
          -69},{72,-24.5},{59.7,-24.5}}, color={0,0,127}));
  connect(senVel.v, eff.v_a) annotation (Line(points={{54,-69},{54,-64},{70,-64},
          {70,-14},{53.5,-14},{53.5,-19.9}}, color={0,0,127}));
  connect(senWetBul.T, eff.Ta_wb1) annotation (Line(points={{28,-69},{28,-60},{
          68,-60},{68,-27.9},{59.7,-27.9}}, color={0,0,127}));
  connect(senTem_a2.T, eff.Ta_db1) annotation (Line(points={{87,-50},{96,-50},{
          96,-31.3},{59.7,-31.3}}, color={0,0,127}));
  connect(senTem_a1.T, eff.Tw1) annotation (Line(points={{-81,90},{-76,90},{-76,
          100},{78,100},{78,-34.9},{59.7,-34.9}}, color={0,0,127}));
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
