within EvaporativeCooler.BaseClasses;
record Breezair_icon_210_RPM560 "System performance curves for Breezair icon 210 under speed 9 (RPM 560)"
  extends Buildings.Fluid.Movers.Data.Generic(
    use_powerCharacteristic=true,
    speed_rpm_nominal=560,
    power(V_flow={1.53936,1.60388,1.68131,1.77380,1.88134,1.98242,2.07276,2.12223}, P={1059.89,1069.94,1090.034,1112.14,1140.27,1162.39,1186.50,1200.56}),
    pressure(V_flow={1.534,1.597,1.653,1.730,1.808,1.874,1.944,2.021,2.092,2.126}, dp={174.188,159.2576,144.3272,121.9316,99.536,79.6288,57.2332,34.8376,12.422}));
end Breezair_icon_210_RPM560;
