function Norm, Axis, Frac
  Coord = Axis.crange[0] * (1d - Frac) + Axis.crange[1] * Frac
  return, (Axis.type EQ 1) ? 10d^Coord : Coord
end
