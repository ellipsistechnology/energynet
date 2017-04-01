package ellipsis.energy.test;

import static ellipsis.util.ArrayHelper.array;
import static ellipsis.util.VectorHelper.vector;
import static ellipsis.util.VectorHelper.vectorArray;

import java.util.Arrays;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.util.GridlessTransformerRegulator;

public class GridlessTrasformerRegulatorTest extends GridTestHelper
{
	private static final int MIN_STEP = -10;
	private static final int MAX_STEP = 10;
	private static final double STEP_SIZE = 0.01;
	private static final int TAP_CHANGE_DELAY = 10;

	public static void main(String[] args)
	{
		assertAssertsOn();
		
		GridlessTrasformerRegulatorTest test = new GridlessTrasformerRegulatorTest();
		test.parseArgs(args);
		
		test.test_schedule();
		test.test_timing();
		test.test_varyingVoltages();
//		test.test_voltageComparison();
	}

//	private void test_voltageComparison()
//	{
//		GridlessTransformerRegulator xfmr = setup();
//		// TODO check voltage
//		double v1 = ;
//		int T = 12;
//		xfmr.setRatioSchedule(array(T, 0.9));
//		double v2 = ;
//		// TODO check voltage has dropped
//		xfmr.setRatioSchedule(array(T, 1.1));
//		double v3 = ;
//		// TODO check vaoltage has risen
//	}

	private GridlessTransformerRegulator setup()
	{
		GridlessTransformerRegulator xfmr = new GridlessTransformerRegulator();
//		xfmr.setBusCount(3);
		xfmr.setBusIndeces(Arrays.asList(1,2,3));
		xfmr.setMinStep(MIN_STEP);
		xfmr.setMaxStep(MAX_STEP);
		xfmr.setStepSize(STEP_SIZE);
		xfmr.setVoltageReference(1);
		xfmr.setTapChangeDelay(TAP_CHANGE_DELAY);
		int T = 12;
		xfmr.setT(T);
//		xfmr.setVoltages(voltageSchedule());
		xfmr.setV0(voltageSchedule()[0]);
		xfmr.setDeltaVExternal(deltaVSchedule());
		xfmr.setRatioSchedule(array(T, 1.0));
		double baseImpedance = 17.3056;
		xfmr.setAdmittance(new Complex(0.001, 0.008).reciprocal().multiply(baseImpedance));
		return xfmr;
	}

	private RealVector[] voltageSchedule()
	{
		return vectorArray(12, 6,
				0, 0, 0, 1.05, 1.05, 1.05,
				0, 0, 0, 1.05, 1.05, 1.05,
				0, 0, 0, 1.06, 1.07, 1.08,
				0, 0, 0, 1.08, 1.07, 1.06,
				0, 0, 0, 1.05, 1.04, 1.03,
				0, 0, 0, 1.0,  1.0,  1.0,
				0, 0, 0, 0.97, 0.95, 0.93,
				0, 0, 0, 0.95, 0.94, 0.93,
				0, 0, 0, 0.92, 0.91, 0.90,
				0, 0, 0, 0.90, 0.91, 0.92,
				0, 0, 0, 0.93, 0.94, 0.95,
				0, 0, 0, 1.0,  1.0,  1.0
		);
	}

	private RealVector[] deltaVSchedule()
	{
		RealVector[] voltageSchedule = voltageSchedule();
		RealVector v0 = voltageSchedule[0];
		RealVector[] deltaV = new RealVector[voltageSchedule.length];
		for (int t = 0; t < deltaV.length; t++)
		{
			deltaV[t] = voltageSchedule[t].subtract(v0);
		}
		return deltaV;
	}

	private double[] correctSchedule()
	{
		double[] correctSchedule = new double[]{
				0.95,
				0.95,
				0.93,
				0.93,
				0.96,
				1.0,
				1.05,
				1.06,
				1.09,
				1.09,
				1.06,
				1.0	
		};
		return correctSchedule;
	}

	private void test_schedule()
	{
		startTest();
		
		// Setup:
		GridlessTransformerRegulator xfmr = setup();
		
		// Schedule:
		xfmr.schedule(null);
		double[] schedule = xfmr.mSchedule;
		double[] correctSchedule = correctSchedule();
		for(int i = 0; i < correctSchedule.length; ++i)
			assertTrue(equals(schedule[i], correctSchedule[i], 1e-12), "Incorrect schedule at time "+i+": "+schedule[i]+", should be "+correctSchedule[i]);
		
		endTest("Schedule");
	}

	private void test_timing()
	{
		startTest();
		
		// Setup:
		GridlessTransformerRegulator xfmr = setup();

		// Schedule:
		xfmr.schedule(null);
		double[] schedule = xfmr.mSchedule;
//		xfmr.applySchedule(schedule);
		
		// Make sure the transformer hasn't stepped yet:
		assertTrue(xfmr.getRatioSchedule()[0] == 1.0, "Tap changed before step delay (no time has passed yet).");
		
		// Pass time, but not enough to step:
		xfmr.passTime(TAP_CHANGE_DELAY/2);
		assertTrue(xfmr.getRatioSchedule()[0] == 1.0, "Tap changed before step delay (only a short time has passed).");
		
		// Pass enough time to step:
		xfmr.passTime(TAP_CHANGE_DELAY/2);
		assertTrue(xfmr.getRatioSchedule()[0] != 1.0, "Tap has not changed after time has passed.");
		
		// Make sure it's only taken one step:
		assertTrue(equals(xfmr.getRatioSchedule()[0], 1.0-STEP_SIZE, 1e-12), "Tap did not take a single step or did not step in the correct direction or did not step the correct amount.");
		
		// Pass more time and make sure the rest of the steps are taken:
		double m = xfmr.getRatioSchedule()[0];
		double target = schedule[0];
		while(m > target)
		{
			assertTrue(equals(xfmr.getRatioSchedule()[0], m, 1e-12), "Tap did not take a single step or did not step in the correct direction or did not step the correct amount.");
			xfmr.passTime(TAP_CHANGE_DELAY);
			m -= STEP_SIZE;
		}
		
		// Make sure it's stepped to the correct position:
		assertTrue(equals(target, m, 1e-12), "Tap changed passed target: m="+m+", target="+target+".");
		
		endTest("Timing");
	}

	private void test_varyingVoltages()
	{
		startTest();
		
		// Setup:
		GridlessTransformerRegulator xfmr = setup();

		// Schedule:
		xfmr.schedule(null);
		double[] schedule = xfmr.mSchedule;
		xfmr.applySchedule(schedule);
		
		// Make sure the transformer hasn't stepped yet:
		assertTrue(xfmr.getRatioSchedule()[0] == 1.0, "Tap changed before step delay (no time has passed yet).");
		
		// Pass time, but not enough to step:
		xfmr.passTime(TAP_CHANGE_DELAY/2);
		assertTrue(xfmr.getRatioSchedule()[0] == 1.0, "Tap changed before step delay (only a short time has passed).");
		
		// Change voltages and reschedule:
		RealVector[] vs = voltageSchedule();
		vs[0] = vector(0, 0, 0, 1.0, 1.0, 1.0); // Should set the tap back to 1.0 after appropriate delay.
//		xfmr.setVoltages(vs);
		xfmr.setV0(vs[0]);
		xfmr.schedule(null);
		double[] newSchedule = xfmr.mSchedule;
		xfmr.applySchedule(newSchedule);
		
		// Make sure it's stepped according to the new schedule:
		xfmr.passTime(TAP_CHANGE_DELAY/2);
		assertTrue(xfmr.getRatioSchedule()[0] == 1.0, "Tap didn't change back to 1.0 after voltage change.");
		
		endTest("Varying voltages");
	}
}