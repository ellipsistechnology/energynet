package ellipsis.energy.grid.res;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.grid.DistributedSource;
import ellipsis.energy.smartgrid.global_local.Forecast;

public abstract class RenewableGenerator extends DistributedSource
{
    public abstract Complex getExpectedPower(Forecast forecast);
    public abstract Complex getExpectedMinimumPower(Forecast forecast);
    public abstract Complex getExpectedMaximumPower(Forecast forecast);
}
