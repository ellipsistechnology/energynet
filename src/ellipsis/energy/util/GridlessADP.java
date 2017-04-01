package ellipsis.energy.util;

import org.apache.commons.math3.linear.RealVector;

public interface GridlessADP extends ADPInterface
{
	GridlessData getGridlessData();
	int stateDimension();
	int voltageDimension();
	int powerDimension();
	int controlDimension();
	int unitCount();
	RealVector powerFromControlAndNoise(RealVector u, RealVector w);
	RealVector voltagesFromControlAndNoise(int t, RealVector u, RealVector w);
	RealVector voltageFromPower(int t, RealVector power);

	int dgPowerOffset();
	int storagePowerOffset();
	int loadPowerOffset();
}