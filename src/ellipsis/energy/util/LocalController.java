package ellipsis.energy.util;

import java.util.List;

import org.apache.commons.math3.linear.RealVector;

/**
 * A generic local controller for use by {@link CentralCoordinator}
 * @author bmillar
 *
 */
public interface LocalController
{
	RealVector[] schedule(double[] ctgs);
	
	/**
	 * The terminating time (horizon).
	 * @return
	 */
	int getT();
	
	void train(double[][] defaultBandwidth);
	
	/**
	 * Mean control (typically zero).
	 * @return
	 */
	RealVector uZero();

	int unitCount();
	
	// TODO possibly move the following to a gridless version
	void setPower0(int k, double p);
	void setDeltaVExternal(RealVector[] deltaVExt);
	void setDeltaVExternal(int t, int i, double deltaV);
	RealVector[] getDeltaVExternal();
	void setV0(RealVector arrayRealVector);
	void setV0(int i, double abs);
	RealVector getV0();
	void setPower0(RealVector arrayRealVector);
	List<Integer> getBusIndeces();
	void setBusIndeces(List<Integer> B);
}