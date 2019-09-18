within EvaporativeCooler.BaseClasses;
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
