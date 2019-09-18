within EvaporativeCooler;
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
    annotation (Placement(transformation(extent={{6,-66},{-6,-54}})));
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
 Received tau1 = " + String(tau1) + "\n");
  assert((massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) or
          tau1 > Modelica.Constants.eps,
"The parameter tau1, or the volume of the model from which tau may be derived, is unreasonably small.
 You need to set massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState to model steady-state.
 Received tau1 = " + String(tau1) + "\n");

 // Check for tau2
  assert((energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) or
          tau2 > Modelica.Constants.eps,
"The parameter tau2, or the volume of the model from which tau may be derived, is unreasonably small.
 You need to set energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState to model steady-state.
 Received tau2 = " + String(tau2) + "\n");
  assert((massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) or
          tau2 > Modelica.Constants.eps,
"The parameter tau2, or the volume of the model from which tau may be derived, is unreasonably small.
 You need to set massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState to model steady-state.
 Received tau2 = " + String(tau2) + "\n");

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
    annotation (Line(points={{28,-60},{6,-60}}, color={0,127,255}));
  connect(senMasFra2.X, water_consumption.W_in) annotation (Line(points={{0,-53.4},
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
  connect(senMasFra2.X, heatTransfer.Wx_1) annotation (Line(points={{0,-53.4},{
          10,-53.4},{10,-52},{22,-52},{22,-24.9667},{13.6,-24.9667}},
                                                                color={0,0,127}));
  connect(senTem2.T, heatTransfer.T_db1) annotation (Line(points={{86,-51.2},{
          86,-52},{82,-52},{82,-21.6667},{13.6,-21.6667}},
                                                  color={0,0,127}));
  connect(preDro2.port_a, Volume_2_air.ports[1])
    annotation (Line(points={{-74,-60},{-48,-60}}, color={0,127,255}));
  connect(Volume_2_air.ports[2], senMasFra2.port_b)
    annotation (Line(points={{-52,-60},{-6,-60}}, color={0,127,255}));
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
</html>", revisions="<html>
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
