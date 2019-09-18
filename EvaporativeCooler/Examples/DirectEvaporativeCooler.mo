within EvaporativeCooler.Examples;
model DirectEvaporativeCooler
  extends Modelica.Icons.Example;
  replaceable package Medium1 = Buildings.Media.Water;
  replaceable package Medium2 = Buildings.Media.Air;

  parameter Modelica.SIunits.MassFlowRate m1_flow_nominal=0.32
    "Nominal mass flow rate";
  parameter Modelica.SIunits.MassFlowRate m2_flow_nominal=3.2
    "Nominal mass flow rate";

  Buildings.Fluid.Sources.Boundary_pT bou_air_indoor(
    redeclare package Medium = Medium2,
    use_T_in=false,
    nPorts=1)
    annotation (Placement(transformation(extent={{-116,-36},{-96,-16}})));
  Modelica.Blocks.Sources.RealExpression Temperature(y=291.70)
    annotation (Placement(transformation(extent={{-142,18},{-122,38}})));
  Modelica.Blocks.Sources.RealExpression watePumpHead(y=0.32)  "in pascals"
    annotation (Placement(transformation(extent={{-142,42},{-122,62}})));
  Buildings.Fluid.Movers.FlowControlled_m_flow Pump(
    redeclare package Medium = Medium1,
    addPowerToMedium=false,
    nominalValuesDefineDefaultPressureCurve=true,
    dp_nominal=11652.14,
    m_flow_start=0.32,
    m_flow_nominal=m1_flow_nominal)
    annotation (Placement(transformation(extent={{-54,12},{-34,32}})));
  Buildings.Fluid.Sources.Boundary_ph Bou_water1(redeclare package Medium =
        Medium1, nPorts=1)
    annotation (Placement(transformation(extent={{80,12},{60,32}})));
  Modelica.Blocks.Sources.RealExpression mflow(y=10*m2_flow_nominal)
    annotation (Placement(transformation(extent={{174,6},{154,-14}})));
  Modelica.Blocks.Sources.RealExpression Air_composition(y=41.2 + 273.15)
    annotation (Placement(transformation(extent={{174,-14},{154,-34}})));
  Buildings.Fluid.Sensors.TemperatureWetBulbTwoPort senWetBul(
    redeclare package Medium = Medium2, m_flow_nominal=m2_flow_nominal)
           annotation (Placement(transformation(extent={{-66,-24},{-56,-14}})));
  Buildings.Fluid.Sensors.MassFlowRate senMasFlo(redeclare package Medium =
        Medium2)
    annotation (Placement(transformation(extent={{-42,-24},{-32,-14}})));
  Buildings.Fluid.Sensors.TemperatureTwoPort senTem1(
    redeclare package Medium = Medium2, m_flow_nominal=m2_flow_nominal)
           annotation (Placement(transformation(extent={{26,-48},{36,-36}})));
  Buildings.Fluid.Sensors.TemperatureWetBulbTwoPort senWetBul1(
    redeclare package Medium = Medium2, m_flow_nominal=m2_flow_nominal)
           annotation (Placement(transformation(extent={{42,-48},{54,-36}})));
  Buildings.Fluid.Sensors.RelativeHumidityTwoPort senRelHum(
    redeclare package Medium = Medium2, m_flow_nominal=m2_flow_nominal)
           annotation (Placement(transformation(extent={{-24,-24},{-14,-14}})));
  Modelica.Blocks.Sources.Ramp ramp1(
    height=2.54,
    duration=10000,
    offset=0.201,
    startTime=100)
    annotation (Placement(transformation(extent={{140,56},{120,76}})));
  Buildings.Fluid.Sensors.TemperatureTwoPort senTem2(
    redeclare package Medium = Medium2, m_flow_nominal=m2_flow_nominal)
           annotation (Placement(transformation(extent={{-80,-24},{-70,-12}})));
  Buildings.Fluid.Sources.Boundary_pT bou_air_outdoor(
    redeclare package Medium = Medium2,
    use_Xi_in=true,
    use_T_in=true,
    nPorts=1)
    annotation (Placement(transformation(extent={{132,-60},{112,-40}})));
  Buildings.Fluid.Sources.Boundary_pT Bou_water2(
    redeclare package Medium = Medium1,
    use_T_in=true,
    nPorts=1)
    annotation (Placement(transformation(extent={{-104,12},{-84,32}})));
  EvaporativeCooler.DirectEvaporativeCooler evaCoo(
    redeclare package Medium1 = Medium1,
    redeclare package Medium2 = Medium2,
    allowFlowReversal1=true,
    allowFlowReversal2=true,
    dp1_nominal=11652.14,
    dp2_nominal=50,
    Thickness=0.1,
    Length=0.50,
    Height=0.50,
    K_value=0.04,
    Contact_surface_area=275,
    Rcon=10,
    m1_flow_nominal=m1_flow_nominal,
    m2_flow_nominal=m2_flow_nominal)
    annotation (Placement(transformation(extent={{2,6},{22,26}})));
  Modelica.Blocks.Tables.CombiTable1D combiTable1D1(table=[0,0.0286; 100,0.0227;
        200,0.0175; 300,0.0129; 400,0.00882; 500,0.00519; 600,0.00193; 700,
        0.0405; 800,0.0329; 900,0.0262; 1000,0.0203; 1100,0.0151; 1200,0.0106;
        1300,0.0065; 1400,0.0198; 1500,0.0152; 1600,0.011; 1700,0.0075; 1800,
        0.00419; 1900,0.00127; 2000,0.00127])
    annotation (Placement(transformation(extent={{204,-102},{184,-82}})));
  Modelica.Blocks.Sources.Clock clock
    annotation (Placement(transformation(extent={{228,-18},{208,2}})));
  Modelica.Blocks.Sources.RealExpression Air_composition1(y=0.0150)
    annotation (Placement(transformation(extent={{178,-44},{158,-64}})));
  Modelica.Blocks.Sources.Ramp ramp2(
    height=2.5,
    duration=10000,
    offset=0,
    startTime=100)
    annotation (Placement(transformation(extent={{208,60},{188,80}})));
  Buildings.Fluid.Movers.FlowControlled_m_flow Fan(
    redeclare package Medium = Medium2,
    y_start=1,
    per(
      pressure(V_flow={1.63,1.67,1.71,1.74,1.78,1.83,1.90,1.95,2.003,2.04,2.09,
            2.14,2.18,2.23}, dp={145.9,136.04,124.4,120,110.1,99.0,85.4,75.55,
            65.67,54.56,42.22,28.64,18.76,7.65}),
      use_powerCharacteristic=true,
      power(V_flow={0.201,2.548516}, P={70,1870}),
      motorCooledByFluid=true,
      speed_rpm_nominal=420),
    inputType=Buildings.Fluid.Types.InputType.Continuous,
    m_flow_start=0.201,
    m_flow_nominal=m2_flow_nominal,
    addPowerToMedium=false)
    annotation (Placement(transformation(extent={{88,-52},{68,-32}})));

equation
  connect(senWetBul.port_b,senMasFlo. port_a)
    annotation (Line(points={{-56,-19},{-42,-19}}, color={0,127,255}));
  connect(senTem1.port_b,senWetBul1. port_a)
    annotation (Line(points={{36,-42},{42,-42}}, color={0,127,255}));
  connect(senMasFlo.port_b,senRelHum. port_a)
    annotation (Line(points={{-32,-19},{-24,-19}}, color={0,127,255}));
  connect(senWetBul.port_a, senTem2.port_b) annotation (Line(points={{-66,-19},
          {-68,-19},{-68,-18},{-70,-18}}, color={0,127,255}));
  connect(senTem2.port_a, bou_air_indoor.ports[1]) annotation (Line(points={{-80,-18},
          {-86,-18},{-86,-26},{-96,-26}},          color={0,127,255}));
  connect(Temperature.y, Bou_water2.T_in) annotation (Line(points={{-121,28},{
          -114,28},{-114,26},{-106,26}}, color={0,0,127}));
  connect(evaCoo.port_b1, Bou_water1.ports[1])
    annotation (Line(points={{22,22},{60,22}}, color={0,127,255}));
  connect(evaCoo.port_a2, senTem1.port_a) annotation (Line(points={{22,10},{24,10},
          {24,-42},{26,-42}}, color={0,127,255}));
  connect(Bou_water2.ports[1], Pump.port_a)
    annotation (Line(points={{-84,22},{-54,22}}, color={0,127,255}));
  connect(Pump.port_b, evaCoo.port_a1) annotation (Line(points={{-34,22},{-16,22},
          {-16,22},{2,22}}, color={0,127,255}));
  connect(watePumpHead.y, Pump.m_flow_in) annotation (Line(points={{-121,52},{
          -88,52},{-88,54},{-44,54},{-44,34}}, color={0,0,127}));
  connect(senRelHum.port_b, evaCoo.port_b2) annotation (Line(points={{-14,-19},{
          -6,-19},{-6,10},{2,10}}, color={0,127,255}));
  connect(clock.y, combiTable1D1.u[1]) annotation (Line(points={{207,-8},{198,
          -8},{198,-92},{206,-92}}, color={0,0,127}));
  connect(Air_composition.y, bou_air_outdoor.T_in)
    annotation (Line(points={{153,-24},{153,-46},{134,-46}}, color={0,0,127}));
  connect(bou_air_outdoor.Xi_in[1], Air_composition1.y) annotation (Line(points=
         {{134,-54},{146,-54},{146,-54},{157,-54}}, color={0,0,127}));
  connect(senWetBul1.port_b, Fan.port_b)
    annotation (Line(points={{54,-42},{68,-42}}, color={0,127,255}));
  connect(Fan.port_a, bou_air_outdoor.ports[1]) annotation (Line(points={{88,
          -42},{100,-42},{100,-50},{112,-50}}, color={0,127,255}));
  connect(mflow.y, Fan.m_flow_in)
    annotation (Line(points={{153,-4},{78,-4},{78,-30}}, color={0,0,127}));
  annotation (
    Diagram(coordinateSystem(extent={{-160,-100},{240,100}})));
end DirectEvaporativeCooler;
