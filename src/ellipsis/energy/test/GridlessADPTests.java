package ellipsis.energy.test;

import static ellipsis.util.ArrayHelper.array;
import static ellipsis.util.VectorHelper.vector;

import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.util.ADP;
import ellipsis.energy.util.CostFunction;
import ellipsis.energy.util.DiscreteGridlessADP;
import ellipsis.util.NonParametricApproximator;
import ellipsis.util.Pair;
import ellipsis.util.VectorHelper;

public class GridlessADPTests extends GridTestHelper
{
	// See VariousGrids.grid01():
	private static final Array2DRowRealMatrix sensitivities = new Array2DRowRealMatrix(new double[][] {
			{0.031155840775389045,	0.031160889802630452,	0.031172191211908975,	-0.06179983216867504,	-0.061799807748500905,	-0.061794291724901855,	},
			{0.031803052148664676,	0.06438691678423646,	0.03187100876952077,	-0.06147872657548379,	-0.1258092232419216,	-0.061455699065703875,	},
			{0.03183791768583748,	0.03185976122237478,	0.06534363174121541,	-0.06146142848793761,	-0.06146132283927689,	-0.12742533560915956,	},
			{0.06421755470365582,	0.06579875808836139,	0.06933801936558166,	0.031860713281419574,	0.03186836094505047,	0.03359581353500421,	},
			{0.06485330528275499,	0.12972752738469534,	0.07002446228868395,	0.03217613274908146,	0.06344090180795961,	0.03392841040837863,	},
			{0.06661728744893677,	0.06825757850144185,	0.1386193769083688,	0.033051309798248295,	0.03305924324586652,	0.06778374971657394,	},
	});
	
	private static final RealVector v_0 = vector(-0.0031, -0.0079, -0.0081, 0.9693, 0.9599, 0.9356); // Args, Abs's
	
	//                                           SOC
	private static final RealVector x_0 = vector(0.5);
	
	//                                               P_DG P_S    P_L          Q_L
	private static final RealVector power_0 = vector(0.2, -0.15, -0.45, 0, 0, -0.15);
	
	private static final RealVector v_t = vector(-0.0015, -0.0059, -0.0047, 0.9811, 0.9724, 0.9592); // Args, Abs's
	
	//                                           SOC
//	private static final RealVector x_t = vector(0.5);
	
	//                                               P_DG  P_S    P_L         Q_L
	private static final RealVector power_t = vector(0.19, -0.14, -0.3, 0, 0, -0.1);
	
	public static void main(String[] args)
	{
		// Ensure that asserts are switched on:
		assertAssertsOn();
		
		GridlessADPTests tests = new GridlessADPTests();
		
		tests.parseArgs(args);
		
		tests.test_enumerateControls();
		tests.test_randomVariations();
		tests.test_f();
		tests.test_g();
		tests.test_voltages();
		tests.test_argMin();
		tests.test_train(0);
		tests.test_external();
		
		if(tests.iterations > 0)
		{
			double averageError = 0;
			for(int seed = 0; seed < tests.iterations; ++seed)
			{
				System.out.print("Iteration "+seed+": ");
//				GridlessADP.rand = new Random(i); Moved to inside test_train()
				averageError += tests.test_train(seed);
			}
			averageError /= tests.iterations;
			System.out.println("Average error = "+averageError);
		}
	}

	private void test_enumerateControls()
	{
		startTest();
		
		DiscreteGridlessADP adp = new DiscreteGridlessADP();
		adp.setX0(x_0);
		
		// Test with default charge/discharge rates and 1 DG 1 storage:
		Set<RealVector> U = adp.enumerateControls(0, x_0);
		assertTrue( U.size() == 6, "Test with default charge/discharge rates and 1 DG 1 storage: U[t].size() = "+U.size()+", expected 6");

		assertTrue( U.contains(vector(0.0,  0.5)) , "Test with default charge/discharge rates and 1 DG 1 storage: U did not contain (0,  0.5)");
		assertTrue( U.contains(vector(0.0,  0.0)) , "Test with default charge/discharge rates and 1 DG 1 storage: U did not contain (0,  0)");
		assertTrue( U.contains(vector(0.0, -1.0)) , "Test with default charge/discharge rates and 1 DG 1 storage: U did not contain (0,  -1)");
		assertTrue( U.contains(vector(1.0,  0.5)) , "Test with default charge/discharge rates and 1 DG 1 storage: U did not contain (1,  0.5)");
		assertTrue( U.contains(vector(1.0,  0.0)) , "Test with default charge/discharge rates and 1 DG 1 storage: U did not contain (1,  0)");
		assertTrue( U.contains(vector(1.0, -1.0)) , "Test with default charge/discharge rates and 1 DG 1 storage: U did not contain (1,  -1)");
		
		RealVector x_t = x_0;
		for (RealVector u : U)
		{
			double P_DG = u.getEntry(0);
			double P_S = u.getEntry(1);
			assertTrue( P_DG >= 0 , "Test with default charge/discharge rates and 1 DG 1 storage: Control contained negative value for DG");
			assertTrue( P_S <= adp.getGridlessData().storageDischargeRate[0] , "Test with default charge/discharge rates and 1 DG 1 storage: Control contained high value for storage: "+P_S);
			assertTrue( P_S >= adp.getGridlessData().storageChargeRate[0] , "Test with default charge/discharge rates and 1 DG 1 storage: Control contained low value for storage: "+P_S);
			double soc_t = x_0.getEntry(0);
			double soc_next = soc_t - P_S;
			assertTrue( soc_next <= adp.getGridlessData().storageMaxCapacity[0] , "Test with default charge/discharge rates and 1 DG 1 storage: Control will lead to over-capacity for storage: "+P_S+", soc_t = "+x_t.getEntry(0));
			assertTrue( soc_next >= 0 , "Test with default charge/discharge rates and 1 DG 1 storage: Control will lead to under-capacity for storage: "+P_S+", soc_t = "+x_t.getEntry(0));
		}
		
		// Test with specified charge/discharge rate and 3 DG and 2 storage:
		adp.getGridlessData().dgCount = 3;
		adp.getGridlessData().storageCount = 2;
		adp.getGridlessData().meanP_DG = new double[adp.T][];
		adp.getGridlessData().storageChargeRate = new double[] {-1.0, -1.5};
		adp.getGridlessData().storageDischargeRate = new double[] {1.0, 1.5};
		adp.getGridlessData().storageMaxCapacity = new double[] {2.0, 3.0};
		for(int t = 0; t < adp.getGridlessData().meanP_DG.length; ++t)
		{
			adp.getGridlessData().meanP_DG[t] = new double[]{1.0, 1.5, 2.0};
		}
		
		RealVector x = vector(0.5, 1); // Only partial state since only SOC components are relevant.
		U = adp.enumerateControls(0, x);
		assertTrue( U.size() == 72 , "Test with specified charge/discharge rate and 3 DG and 2 storage: U[t].size() = "+U.size()+", expected 72");
		
		assertTrue( U.contains(vector(0.0, 0.0, 0.0,  0.5,  1.0)) , "Test with specified charge/discharge rate and 3 DG and 2 storage: U did not contain (0.0, 0.0, 0.0, -0.5, -1.0)");
		assertTrue( U.contains(vector(1.0, 0.0, 0.0, -1.0, -1.5)) , "Test with specified charge/discharge rate and 3 DG and 2 storage: U did not contain (1.0, 0.0, 0.0,  1.0,  1.5)");
		assertTrue( U.contains(vector(1.0, 1.5, 0.0,  0.5,  1.0)) , "Test with specified charge/discharge rate and 3 DG and 2 storage: U did not contain (1.0, 1.5, 0.0, -0.5, -1.0)");
		assertTrue( U.contains(vector(0.0, 1.5, 0.0, -1.0,  0.0)) , "Test with specified charge/discharge rate and 3 DG and 2 storage: U did not contain (0.0, 1.5, 0.0,  1.0,  0.0)");
		assertTrue( U.contains(vector(0.0, 0.0, 2.0,  0.5, -1.5)) , "Test with specified charge/discharge rate and 3 DG and 2 storage: U did not contain (0.0, 0.0, 2.0, -0.5,  1.5)");
		
		endTest("test_enumerateControls()");
	}

	private void test_randomVariations()
	{
		startTest();
		
		DiscreteGridlessADP adp = new DiscreteGridlessADP();
//		RealVector u_0 = vector(1.0, 0.0); // {P_DG, P_S}
		
		RealVector w_0a = adp.randomVariations(0);//, u_0);
		assertTrue( w_0a != null , "randomVariations() returned null");
		
		RealVector w_0b = adp.randomVariations(0);//, u_0);
		assertTrue( w_0b != null , "randomVariations() returned null");
		assertTrue( w_0a.getDimension() == w_0b.getDimension() , "randomVariations() successive calls returned vectors of different lengths");
		assertTrue( !w_0a.equals(w_0b) , "randomVariations() successive calls returned the same vector");
		
		// Test random variations mean and S.D:
		RealVector meanSum = new ArrayRealVector(w_0a.getDimension());
		int samples = 1000000;
		RealVector[] w = new RealVector[samples];
		for(int i = 0; i < samples; ++i)
		{
			w[i] = adp.randomVariations(0);//, u_0);
			meanSum = meanSum.add(w[i]);
		}
		RealVector mean = meanSum.mapDivide(samples);
		assertTrue( equals(mean.getEntry(0), 0, 1e-3) , "Test random variations mean and S.D: Mean of DG variations was not 0, was "+mean.getEntry(0));
		assertTrue( equals(mean.getEntry(1), adp.getGridlessData().meanP_L[0][0], 1e-3) , "Test random variations mean and S.D: Mean of P_L variations was not "+adp.getGridlessData().meanP_L[0][0]+", was "+mean.getEntry(1));
		assertTrue( equals(mean.getEntry(2), adp.getGridlessData().meanQ_L[0][0], 1e-3) , "Test random variations mean and S.D: Mean of Q_L variations was not "+adp.getGridlessData().meanQ_L[0][0]+", was "+mean.getEntry(2));
		
		RealVector varSum = new ArrayRealVector(w_0a.getDimension());
		for(int i = 0; i < samples; ++i)
		{
			RealVector sd = w[i].subtract(mean);
			RealVector var = sd.ebeMultiply(sd);
			varSum = varSum.add(var);
		}
		RealVector var = varSum.mapDivide(samples);
		RealVector sd = new ArrayRealVector(var.getDimension());
		for(int i = 0; i < var.getDimension(); ++i)
		{
			double d = var.getEntry(i);
			sd.setEntry(i, Math.sqrt(d));
		}
		assertTrue( equals(sd.getEntry(0), adp.getGridlessData().sdP_DG[0][0], 1e-2) , "Test random variations mean and S.D: S.D. of DG variations was not "+adp.getGridlessData().sdP_DG[0][0]);
		assertTrue( equals(sd.getEntry(1), adp.getGridlessData().sdP_L[0][0], 1e-2) , "Test random variations mean and S.D: S.D. of P_L variations was not "+adp.getGridlessData().sdP_L[0][0]);
		assertTrue( equals(sd.getEntry(2), adp.getGridlessData().sdQ_L[0][0], 1e-2) , "Test random variations mean and S.D: S.D. of Q_L variations was not "+adp.getGridlessData().sdQ_L[0][0]);
		
		endTest("test_randomVariations()");
	}

	private void test_f()
	{
		startTest();
		
		DiscreteGridlessADP adp = new DiscreteGridlessADP();
		adp.setX0(x_0);
		
		// Test with no control and no noise:
		RealVector u_noNoise = power_0.getSubVector(0, adp.getGridlessData().dgCount+adp.getGridlessData().storageCount); // power_0 is [P_DG P_S P_L Q_DG Q_S Q_L], u is [P_DG P_S]
		RealVector w_noNoise = adp.wZero(0);
		RealVector x_1 = adp.f(0, x_0, u_noNoise, w_noNoise);
		
		assertTrue( x_1.getEntry(0) == x_0.getEntry(0)-u_noNoise.getEntry(1) , "Test with no control and no noise: x_1[0] != x_0[0] + u_0[1] = 0.5 + 0.15 = 0.65");
		for(int i = 1; i < x_0.getDimension(); ++i)
		{
			assertTrue( x_1.getEntry(i) == x_0.getEntry(i) , "Test with no control and no noise: x_1["+i+"] != x_0["+i+"]; "+x_1.getEntry(i)+" != "+x_0.getEntry(i));
		}
			
		// Test with noise but no control:
		// No longer relevant since state is independent of noise.
//		RealVector w_noise = vector(-0.1, x_1.getEntry(3)+0.1, x_1.getEntry(6)-0.1);
//		x_1 = adp.f(x_0, u_noNoise, w_noise);
//		int i = 0;
//		assertTrue( x_1.getEntry(0) == x_0.getEntry(0)-u_noNoise.getEntry(1) , "Test with noise but no control: x_1[0] != x_0[0] + u_0[1] = 0.5 + 0.15 = 0.65"); // SOC
//		++i;
//		assertTrue( x_1.getEntry(i) == x_0.getEntry(i)-0.1 , "Test with noise but no control: x_1["+i+"] != x_0["+i+"]-0.1"); // P_DG
//		++i;
//		assertTrue( x_1.getEntry(i) == x_0.getEntry(i) , "Test with noise but no control: x_1["+i+"] != x_0["+i+"]"); // P_S
//		++i;
//		assertTrue( x_1.getEntry(i) == x_0.getEntry(i)+0.1 , "Test with noise but no control: x_1["+i+"] != x_0["+i+"]+0.1"); // P_L
//		++i;
//		assertTrue( x_1.getEntry(i) == 0.0 , "Test with noise but no control: x_1["+i+"] != 0.0"); // Q_DG
//		++i;
//		assertTrue( x_1.getEntry(i) == 0.0 , "Test with noise but no control: x_1["+i+"] != 0.0"); // Q_S
//		++i;
//		assertTrue( x_1.getEntry(i) == x_0.getEntry(i)-0.1 , "Test with noise but no control: x_1["+i+"] != x_0["+i+"]-0.1"); // Q_L
		
		// Test with control:
		RealVector u = vector(1.0, 0.5); // P_DG P_S
		x_1 = adp.f(0, x_0, u, w_noNoise);
		int i = 0;
		assertTrue( x_1.getEntry(i) == x_0.getEntry(i)-0.5 , "Test with control: x_1["+i+"] != x_0["+i+"] + 0.5"); // SOC
//		++i;
//		assertTrue( equals(x_1.getEntry(i), x_0.getEntry(i), 1e-6) , "Test with noise and control: x_1["+i+"] != x_0["+i+"]; "+x_1.getEntry(i)+" != "+x_0.getEntry(i)); // P_DG
//		++i;
//		assertTrue( x_1.getEntry(i) == 0.5 , "Test with noise and control: x_1["+i+"] != 0.5"); // P_S
//		++i;
//		assertTrue( x_1.getEntry(i) == x_0.getEntry(i)+0.1 , "Test with noise and control: x_1["+i+"] != x_0["+i+"]+0.1"); // P_L
//		++i;
//		assertTrue( x_1.getEntry(i) == 0.0 , "Test with noise and control: x_1["+i+"] != 0.0"); // Q_DG
//		++i;
//		assertTrue( x_1.getEntry(i) == 0.0 , "Test with noise and control: x_1["+i+"] != 0.0"); // Q_S
//		++i;
//		assertTrue( x_1.getEntry(i) == x_0.getEntry(i)-0.1 , "Test with noise and control: x_1["+i+"] != x_0["+i+"]-0.1"); // Q_L
		
		endTest("test_f()");
	}

	private void test_g()
	{
		startTest();
		
		DiscreteGridlessADP adp = new DiscreteGridlessADP();
		adp.setCostFunction(new VoltageCostFunction());
		adp.setX0(x_0);
		adp.getGridlessData().power_0 = power_0;
		RealVector v_0 = vector(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
		adp.getGridlessData().v_0 = v_0; // Complex values for DG, S and L busses (args and abs's respectively).
		adp.getGridlessData().sensitivities = sensitivities;
		adp.getGridlessData().meanP_L = new double[][]{{power_0.getEntry(2)}, {0}};
		adp.getGridlessData().meanQ_L = new double[][]{{power_0.getEntry(5)}, {0}};
		adp.getGridlessData().setDeltaVExternal(vector(6, 0.0));
		
		// Test voltage cost:
		double vc_a = voltageCost(1);
		assertTrue( vc_a == 0 , "Test voltage cost: voltage cost of 1 != 0, was "+vc_a);
		double vc_b = voltageCost(1.01);
		assertTrue( vc_b > 0 , "Test voltage cost: voltage cost of 1.01 <= 0, was "+vc_b);
		double vc_c = voltageCost(1.05);
		assertTrue( 1-vc_c < 1e-3 , "Test voltage cost: voltage cost of 1.05 != 1, was "+vc_c);
		double vc_d = voltageCost(0.95);
		assertTrue( 1-vc_d < 1e-3 , "Test voltage cost: voltage cost of 0.95 != 1, was "+vc_d);
		
		// Test at x_0 with no controls:
		RealVector u_noNoise = vector(power_0.getEntry(0), power_0.getEntry(1));
		double g = adp.g(0, x_0, u_noNoise);
		double cost = 0;
		for(int i = v_0.getDimension()/2; i < v_0.getDimension(); ++i)
		{
			cost += voltageCost(v_0.getEntry(i));
		}
		cost /= v_t.getDimension()/2;
		assertTrue( g == cost , "Test at x_0 with no controls: cost is not as expected, should be cost of voltage in state x_0.");
		
		// Test x_t with controls:
		RealVector x_t = vector(0.6);
		RealVector power_t = vector(0.4, 0.0, -0.5, 0.0, 0.0, -0.1);
		adp.getGridlessData().meanP_L[1] = new double[] {power_t.getEntry(2)};
		adp.getGridlessData().meanQ_L[1] = new double[] {power_t.getEntry(5)};
//		RealVector deltaX_t = x_t.subtract(x_0);
		RealVector deltaS_t = power_t.subtract(power_0);
//		RealVector deltaS_t = deltaX_t.getSubVector(1, deltaX_t.getDimension()-1);
		RealVector deltaV_t = adp.getGridlessData().sensitivities.operate(deltaS_t);
		RealVector v_t = v_0.add(deltaV_t);
		RealVector u_t = vector(power_t.getEntry(0), power_t.getEntry(1));
				
		g = adp.g(1, x_t, u_t); // t=1 due to load means specified above
		
		cost = 0;
		for(int i = v_t.getDimension()/2; i < v_t.getDimension(); ++i)
		{
			cost += voltageCost(v_t.getEntry(i));
		}
		cost /= v_t.getDimension()/2;
		assertTrue( g == cost , "Test x_t with controls: cost is not as expected, should be average cost of voltages in state x_t.");
		
		endTest("test_g()");
	}

	public static double voltageCost(double v)
	{
		double exp = -Math.pow(3*(1 - v)/0.05, 2);
		double cost = 1 - Math.exp(exp);
		return cost;
	}
	
	public static class VoltageCostFunction implements CostFunction
	{
		@Override
		public double g(ADP adp, int t, RealVector x, RealVector u)
		{
			if(u == null)
				return 0;
			
			double sum = 0;
			RealVector v = ((DiscreteGridlessADP)adp).voltagesFromControlAndNoise(t, u, adp.wZero(t));
			int fullDimension = v.getDimension();
			int dimension = fullDimension/2;
			for(int i = dimension; i < fullDimension; ++i)
			{
				double cost = voltageCost(v.getEntry(i));
				sum += cost;
			}
			return sum/dimension;
		}
	}

	private void test_voltages()
	{
		startTest();
		
		DiscreteGridlessADP adp = new DiscreteGridlessADP();
		adp.getGridlessData().power_0 = power_0;
		setupADP(adp);
		
		// Test voltages at x_0:
		RealVector v = adp.voltageFromPower(0, power_0);
		for(int i = 0; i < v_0.getDimension(); ++i)
		{
			assertTrue( v.getEntry(i) == v_0.getEntry(i) , "Test voltages at x_0: v["+i+"] != v_0["+i+"]; v["+i+"] = "+v.getEntry(i));
		}
		
		// Test voltages at x_t:
		v = adp.voltageFromPower(0, power_t);
		for(int i = 0; i < v_t.getDimension()/2; ++i)
		{
			assertTrue( equals(v.getEntry(i), v_t.getEntry(i), 5e-4) , "Test voltages at x_t: v["+i+"] != v_t["+i+"]; v["+i+"] = "+v.getEntry(i)+", v_t["+i+"] = "+v_t.getEntry(i));
		}
		for(int i = v_t.getDimension()/2; i < v_t.getDimension(); ++i)
		{
			assertTrue( equals(v.getEntry(i), v_t.getEntry(i), 5e-3) , "Test voltages at x_t: v["+i+"] != v_t["+i+"]; v["+i+"] = "+v.getEntry(i)+", v_t["+i+"] = "+v_t.getEntry(i));
		}
		
		// Test that increased DG produces higher voltages:
		{
			RealVector power_t2 = new ArrayRealVector(power_t);
			power_t2.setEntry(1, power_t.getEntry(0)+0.1);
			RealVector v2 = adp.voltageFromPower(0, power_t2);
			for(int i = 0; i < v2.getDimension(); ++i)
			{
				assertTrue( v2.getEntry(i) > v.getEntry(i) , "Test that increased DG produces higher voltages: At bus "+i+", v2="+v2.getEntry(i)+" <= "+v.getEntry(i));
			}
		}
		
		// Test that decreased DG produces lower voltages:
		{
			RealVector power_t2 = new ArrayRealVector(power_t);
			power_t2.setEntry(0, power_t.getEntry(0)-0.1);
			RealVector v2 = adp.voltageFromPower(0, power_t2);
			for(int i = 0; i < v2.getDimension(); ++i)
			{
				assertTrue( v2.getEntry(i) < v.getEntry(i) , "Test that decreased DG produces lower voltages: At bus "+i+", v2="+v2.getEntry(i)+" <= "+v.getEntry(i));
			}
		}
		
		// Test that increased storage produces higher voltages:
		{
			RealVector power_t2 = new ArrayRealVector(power_t);
			power_t2.setEntry(2, power_t.getEntry(2)+0.1);
			RealVector v2 = adp.voltageFromPower(0, power_t2);
			for(int i = 0; i < v2.getDimension(); ++i)
			{
				assertTrue( v2.getEntry(i) > v.getEntry(i) , "Test that increased storage produces higher voltages: At bus "+i+", v2="+v2.getEntry(i)+" >= "+v.getEntry(i));
			}
		}
		
		// Test that decreased storage produces lower voltages:
		{
			RealVector power_t2 = new ArrayRealVector(power_t);
			power_t2.setEntry(2, power_t.getEntry(2)-0.1);
			RealVector v2 = adp.voltageFromPower(0, power_t2);
			for(int i = 0; i < v2.getDimension(); ++i)
			{
				assertTrue( v2.getEntry(i) < v.getEntry(i) , "Test that decreased storage produces lower voltages: At bus "+i+", v2="+v2.getEntry(i)+" <= "+v.getEntry(i));
			}
		}
		
		endTest("test_voltages()");
	}

	private void test_argMin()
	{
		startTest();
		
		// Make g() return ||u||^2:
		DiscreteGridlessADP adp = new DiscreteGridlessADP() {
			@Override
			public double g(int t, RealVector x, RealVector u)
			{
				double norm = u.getNorm();
				return norm*norm;
			}
			
			@Override
			public boolean breached(int t, RealVector u, RealVector w)
			{
				return false; // No constraints.
			}
		};
		setupADP(adp);
		
		Set<RealVector> U = adp.enumerateControls(0, x_0);
		RealVector u = adp.argMin(0, x_0, U, adp.wZero(0)); // Should find arg min ||u||^2 = {0.0}.
		for(int i = 0; i < u.getDimension(); ++i)
			assertTrue( u.getEntry(i) == 0.0 , "g() = ||u||^2: u["+i+"] != 0");
		
		// Minimal ||u|| (as above), move SOC from 0.5 to 5:
		{
			adp.approxJ = new NonParametricApproximator[adp.T];
			for(int t = 0; t < adp.T; ++t)
			{
				adp.approxJ[t] = new NonParametricApproximator() {
					@Override
					public double value(RealVector x)
					{
						double diff = x.getEntry(0)-5.0; // SOC deviation from 1.0
						return diff*diff;
					}
				};
			}
			RealVector x = x_0;
			for(int t = 0; t < adp.T; ++t)
			{
				u = adp.argMin(t, x, U, adp.wZero(0)); // Should minimise u and move SOC to 5.0.
				double u_correct = (2*x.getEntry(0)-10)/4; // Analytical solution assuming continuous u.
				assertTrue( u.getEntry(0) == 0 , "ctg = ||u||^2 + ||x-5||^2: u_"+t+"[0] was not 0");
				assertTrue( u.getEntry(1) == Math.max(-1.0, Math.round(u_correct)) , "ctg = ||u||^2 + ||x-5||^2: u_"+t+"[1] = "+u.getEntry(1)+" was not as expected; correct value is "+u_correct);
//				RealVector w_noNoise = vector(0.0, x.getEntry(3), x.getEntry(6));
				x = adp.f(t, x, u, null);//w_noNoise);
			}
		}
		
		// Minimal ||u|| (as above), move SOC from 4.5 to 5:
		{
			RealVector x = new ArrayRealVector(x_0);
			x.setEntry(0, 4.5);
			for(int t = 0; t < adp.T; ++t)
			{
				u = adp.argMin(t, x, U, adp.wZero(0)); // Should minimise u and move SOC to 5.0.
				double u_correct = (2*x.getEntry(0)-10)/4; // Analytical solution assuming continuous u.
				assertTrue( u.getEntry(0) == 0 , "ctg = ||u||^2 + ||x-5||^2: u_"+t+"[0] was not 0");
				assertTrue( u.getEntry(1) == Math.max(-1.0, Math.round(u_correct)) , "ctg = ||u||^2 + ||x-5||^2: u_"+t+"[1] = "+u.getEntry(1)+" was not as expected; correct value is "+u_correct);
//				RealVector w_noNoise = vector(0.0, x.getEntry(3), x.getEntry(6));
				x = adp.f(t, x, u, null);//w_noNoise);
			}
		}
		
		// Make g() return ||u-1||^2:
		adp = new DiscreteGridlessADP() {
			@Override
			public double g(int t, RealVector x, RealVector u)
			{
				RealVector diff = u.subtract(vector(u.getDimension(), 1.0));
				double norm = diff.getNorm();
				return norm*norm;
			}
			
			@Override
			public boolean breached(int t, RealVector u, RealVector w)
			{
				return false; // No constraints.
			}
		};
		setupADP(adp);
		
		RealVector x = new ArrayRealVector(x_0);
		x.setEntry(0, 1);
		U = adp.enumerateControls(0, x);
		u = adp.argMin(0, x, U, adp.wZero(0)); // Should find arg min ||u-1||^2 = {1.0}.
		for(int i = 0; i < u.getDimension(); ++i)
			assertTrue( u.getEntry(i) == 1.0 , "g() return ||u-1||^2: u["+i+"] != 1");
		
		endTest("test_argMin()");
	}

	protected void setupADP(DiscreteGridlessADP adp)
	{
		adp.setX0(x_0);
		adp.getGridlessData().v_0 = v_0; // Complex values for DG, S and L busses (args and abs's respectively).
		adp.getGridlessData().power_0 = power_0;
		adp.getGridlessData().sensitivities = sensitivities;
		adp.getGridlessData().setMeanP_L(array(adp.getT(), -0.5));
		
		// Replace estimators to always return zero:
		adp.approxJ = new NonParametricApproximator[adp.T];
		for(int t = 0; t < adp.T; ++t)
		{
			adp.approxJ[t] = new NonParametricApproximator() {
				@Override
				public double value(RealVector x)
				{
					return 0;
				}
			};
		}
	}
	
	/**
	 * 
	 * @param i 
	 * @return The average error for all times and each test.
	 */
	private double test_train(int seed)
	{
		startTest();
		double localAverageError = 0;
		int errorCount = 0;
		
		// Test training with large storage capacity and zero noise:
		{
			if(verbose)
				System.out.println("\n\nTest training with large storage capacity and zero noise:");
			
			DiscreteGridlessADP adp = setupForTraining(/*withNoise=*/false);
			adp.rand = new Random(seed);
			adp.getGridlessData().storageMaxCapacity = new double[]{10};
			adp.setX0(new ArrayRealVector(adp.getX0()));
			adp.getX0().setEntry(0, 5);
			
			List<Pair<RealVector, Double>> schedule = dp(adp); // Make deterministic DP for comparison:
			
			double[][] bw = new double[adp.T][1];
			for(int t = 0; t < adp.T; ++t)
				bw[t] = new double[]{0.5};
			adp.train(bw);
			double[] approximateCtgs = new double[adp.T+1];
			RealVector[] approximateSchedule = adp.schedule(approximateCtgs);
			
			if(verbose)
			{
				// Approx and true comparison:
				System.out.println("\ncorrect,,,approximate");
				for(int t = 0; t < adp.T; ++t)
				{
					System.out.println(VectorHelper.printVector(schedule.get(t).getKey())+schedule.get(t).getValue()+","+VectorHelper.printVector(approximateSchedule[t])+approximateCtgs[t]);
				}

				// Voltages:
				printVoltageProfiles(adp, schedule);
			}
	
			// Check that it's the same as the deterministic case:
			for(int t = 0; t < adp.T; ++t)
			{
				assertTrue( equals(schedule.get(t).getValue(), approximateCtgs[t], 1) , "Test training with large storage capacity and zero noise: Schedule incorrect at time "+t);
				localAverageError += percentageError(schedule, approximateCtgs, t);
				++errorCount;
			}
		}
		
		// Test training with large storage capacity and with noise:
		{
			if(verbose)
				System.out.println("\n\nTest training with large storage capacity and with noise:");
			
			DiscreteGridlessADP adp = setupForTraining(/*withNoise=*/true);
			adp.getGridlessData().storageMaxCapacity = new double[]{10};
			adp.setX0(new ArrayRealVector(adp.getX0()));
			adp.getX0().setEntry(0, 5);
			
			List<Pair<RealVector, Double>> schedule = dp(adp); // Make deterministic DP for comparison:
			
			double[][] bw = new double[adp.T][1];
			for(int t = 0; t < adp.T; ++t)
				bw[t] = new double[]{0.5};
			adp.train(bw);
			double[] approximateCtgs = new double[adp.T+1];
			RealVector[] approximateSchedule = adp.schedule(approximateCtgs);
			
			if(verbose)
			{
				System.out.println("\ncorrect,,,approximate");
				for(int t = 0; t < adp.T; ++t)
				{
					System.out.println(VectorHelper.printVector(schedule.get(t).getKey())+schedule.get(t).getValue()+","+VectorHelper.printVector(approximateSchedule[t])+approximateCtgs[t]);
				}
			}
	
			for(int t = 0; t < adp.T; ++t)
			{
				assertTrue( equals(schedule.get(t).getValue(), approximateCtgs[t], 1) , "Test training with large storage capacity and with noise: Schedule incorrect at time "+t+"; "+schedule.get(t).getValue()+" != "+approximateCtgs[t]);
				localAverageError += percentageError(schedule, approximateCtgs, t);
				++errorCount;
			}
		}
		
		// Test training with limited storage capacity and with noise:
		{
			if(verbose)
				System.out.println("\n\nTest training with limited storage capacity and with noise:");
			
			DiscreteGridlessADP adp = setupForTraining(/*withNoise=*/true);
			adp.getGridlessData().storageMaxCapacity = new double[]{2}; // Max capacity 2 with initial capacity 0.5.
			
			List<Pair<RealVector, Double>> schedule = dp(adp); // Make deterministic DP for comparison:

			double[][] bw = new double[adp.T][1];
			for(int t = 0; t < adp.T; ++t)
				bw[t] = new double[]{0.5};
			adp.train(bw);
			double[] approximateCtgs = new double[adp.T+1];
			RealVector[] approximateSchedule = adp.schedule(approximateCtgs);
			
			if(verbose)
			{
				System.out.println("\ncorrect,,,approximate");
				for(int t = 0; t < adp.T; ++t)
				{
					System.out.println(VectorHelper.printVector(schedule.get(t).getKey())+schedule.get(t).getValue()+","+VectorHelper.printVector(approximateSchedule[t])+approximateCtgs[t]);
				}
			}
	
			int t = 0; // Only need the cost-to-go at t=0 to be close in this case.
			assertTrue( equals(schedule.get(t).getValue(), approximateCtgs[t], 2) , "Test training with limited storage capacity and with noise: Schedule incorrect at time "+t+", correct CTG is "+schedule.get(t).getValue()+", approximate value is "+approximateCtgs[t]);
			for(t = 0; t < adp.T; ++t)
			{
				localAverageError += percentageError(schedule, approximateCtgs, t);
				++errorCount;
			}
		}
		
		// Test training with limited storage capacity and varying mean DG output:
		{
			if(verbose)
				System.out.println("\n\nTest training with limited storage capacity and varying mean DG output:");
			
			DiscreteGridlessADP adp = setupForTraining(/*withNoise=*/true);
			adp.getGridlessData().storageMaxCapacity = new double[]{2}; // Max capacity 2 with initial capacity 0.5.
			adp.getGridlessData().meanP_DG = new double[][]{{0.2},{0.22},{0.24},{0.22},{0.2},{0.18},{0.15},{0.1},{0.05},{0.0},{0.0},{0.0}}; // [time][unit]
			adp.getGridlessData().sdP_DG = new double[][]{{0.02},{0.02},{0.02},{0.02},{0.02},{0.02},{0.015},{0.015},{0.01},{0.01},{0.0},{0.0}}; // [time][unit]
			
			List<Pair<RealVector, Double>> schedule = dp(adp); // Make deterministic DP for comparison:
			
			double[][] bw = new double[adp.T][1];
			for(int t = 0; t < adp.T; ++t)
				bw[t] = new double[]{0.5};
			adp.train(bw);
			double[] approximateCtgs = new double[adp.T+1];
			RealVector[] approximateSchedule = adp.schedule(approximateCtgs);
			
			if(verbose)
			{
				System.out.println("\ncorrect,,,approximate");
				for(int t = 0; t < adp.T; ++t)
				{
					System.out.println(VectorHelper.printVector(schedule.get(t).getKey())+schedule.get(t).getValue()+","+VectorHelper.printVector(approximateSchedule[t])+approximateCtgs[t]);
				}
			}
	
			int t = 0; // Only need the cost-to-go at t=0 to be close in this case.
			assertTrue( equals(schedule.get(t).getValue(), approximateCtgs[t], 1) , "Test training with limited storage capacity and varying mean DG output: Schedule incorrect at time "+t+", correct CTG is "+schedule.get(t).getValue()+", approximate value is "+approximateCtgs[t]);
			for(t = 0; t < adp.T; ++t)
			{
				localAverageError += percentageError(schedule, approximateCtgs, t);
				++errorCount;
			}
		}
		
		endTest("test_train()");
		
		return localAverageError/errorCount;
	}

	protected DiscreteGridlessADP setupForTraining(final boolean withNoise)
	{
		DiscreteGridlessADP adp = new DiscreteGridlessADP() {
			@Override
			public RealVector[] schedule(double[] ctgs)
			{
				if(verbose)
				{
					int _t = 9;
					System.out.println("\nSamples at t="+_t+":");
					System.out.println("SOC,P_DG,Q_DG,P_S,Q_S,P_L,Q_L,ctg'");
	
					List<Pair<RealVector, Double>> samples = ((NonParametricApproximator)approxJ[_t]).getSamplePoints();
					for (Pair<RealVector, Double> sample : samples)
					{
						RealVector x = sample.getKey();
						double ctg = sample.getValue();
						System.out.println(VectorHelper.printVector(x)+ctg);
					}
	
					System.out.println("\nEstimation at t="+_t+":");
					double dg = 0;
					System.out.println("dg="+dg);
					for(double s = -0.15; s <= 0.15; s += 0.01)
					{
						                    //SOC  P_DG P_S P_L    Q_DG Q_S Q_L
						RealVector x = vector(0.0, dg,  s,  -0.45, 0,   0,  -0.15);
						System.out.println(s+","+approxJ[_t].value(x));
					}
					dg = 0.2;
					System.out.println("dg="+dg);
					for(double s = -0.15; s <= 0.15; s += 0.01)
					{
	                                        //SOC  P_DG P_S P_L    Q_DG Q_S Q_L
						RealVector x = vector(0.0, dg,  s,  -0.45, 0,   0,  -0.15);
						System.out.println(s+","+approxJ[_t].value(x));
					}
				}
				
				return super.schedule(ctgs);
			}
			
			@Override
			public RealVector randomVariations(int t)//, RealVector u_t)
			{
				if(withNoise)
					return super.randomVariations(t);//, u_t);
				else
					return wZero(t);
			}
			
			@Override
			public double g(int t, RealVector x, RealVector u)
			{
				return super.g(t, x, u);
			}
		};
		adp.setCostFunction(new VoltageCostFunction());
		adp.getGridlessData().setDeltaVExternal(vector(6, 0.0));
		
		adp.N_r = 50;
		
		adp.setX0(x_0);
		adp.getGridlessData().power_0 = power_0;
		adp.getGridlessData().v_0 = v_0; // Complex values for DG, S and L busses (args and abs's respectively).
		adp.getGridlessData().sensitivities = sensitivities;
		
		adp.getGridlessData().meanP_L = new double[][]{{-0.4},{-0.4},{-0.4},{-0.4},{-0.4},{-0.4},{-0.4},{-0.4},{-0.4},{-0.4},{-0.4},{-0.4}}; // [time][unit]
		adp.getGridlessData().meanQ_L = new double[][]{{-0.15},{-0.15},{-0.15},{-0.15},{-0.15},{-0.15},{-0.15},{-0.15},{-0.15},{-0.15},{-0.15},{-0.15}}; // [time][unit]
		adp.getGridlessData().sdP_L = new double[][]{{0.03},{0.03},{0.03},{0.03},{0.03},{0.03},{0.03},{0.03},{0.03},{0.03},{0.03},{0.03}}; // [time][unit]
		adp.getGridlessData().sdQ_L = new double[][]{{0.01},{0.01},{0.01},{0.01},{0.01},{0.01},{0.01},{0.01},{0.01},{0.01},{0.01},{0.01}}; // [time][unit]
		
		adp.getGridlessData().storageDischargeRate = new double[]{0.15}; // [unit]
		adp.getGridlessData().storageChargeRate = new double[]{-0.15}; // [unit]
		
		adp.getGridlessData().meanP_DG = new double[][]{{0.2},{0.2},{0.2},{0.2},{0.2},{0.2},{0.2},{0.2},{0.2},{0.2},{0.2},{0.2}}; // [time][unit]
		adp.getGridlessData().sdP_DG = new double[][]{{0.02},{0.02},{0.02},{0.02},{0.02},{0.02},{0.02},{0.02},{0.02},{0.02},{0.02},{0.02}}; // [time][unit]
		
		return adp;
	}

	private void test_external()
	{
		startTest();

		DiscreteGridlessADP adp = new DiscreteGridlessADP();
		adp.getGridlessData().power_0 = power_0;
		setupADP(adp);
		
		// Get comparison with zero external influence:
		RealVector v_base[] = new RealVector[adp.T+1];
		for(int t = 0; t <= adp.T; ++t)
		{
			v_base[t] = adp.voltageFromPower(t, power_t);
		}
		
		// Test \Delta V_ext:
		adp.getGridlessData().deltaV_external = new RealVector[adp.T+1];
		
		RealVector deltaV_ext[] = new RealVector[adp.T+1];
		deltaV_ext[0] = vector(0.01, 0.02, 0.03, 0.4, 0.5, 0.6);
		deltaV_ext[1] = vector(0.02, 0.02, 0.02, 0.4, 0.5, 0.6);
		deltaV_ext[2] = vector(0.03, 0.03, 0.03, 0.4, 0.5, 0.6);
		deltaV_ext[3] = vector(0.01, 0.02, 0.03, 0.4, 0.4, 0.4);
		deltaV_ext[4] = vector(0.01, 0.02, 0.03, 0.5, 0.5, 0.5);
		deltaV_ext[5] = vector(0.01, 0.02, 0.03, 0.6, 0.6, 0.6);
		deltaV_ext[6] = vector(0.02, 0.02, 0.03, 0.4, 0.5, 0.6);
		deltaV_ext[7] = vector(0.02, 0.02, 0.03, 0.5, 0.5, 0.6);
		deltaV_ext[8] = vector(0.01, 0.03, 0.03, 0.4, 0.5, 0.6);
		deltaV_ext[9] = vector(0.01, 0.03, 0.03, 0.5, 0.5, 0.6);
		deltaV_ext[10] = vector(0.01, 0.02, 0.02, 0.4, 0.4, 0.6);
		deltaV_ext[11] = vector(0.01, 0.02, 0.02, 0.4, 0.5, 0.5);
		deltaV_ext[12] = vector(0.01, 0.02, 0.02, 0.4, 0.6, 0.6);
		
		for(int t = 0; t <= adp.T; ++t)
		{
			adp.getGridlessData().deltaV_external[t] = deltaV_ext[t];
		}
		
		RealVector v[] = new RealVector[adp.T+1];
		for(int t = 0; t <= adp.T; ++t)
		{
			v[t] = adp.voltageFromPower(t, power_t);
			
			assertTrue(v[t] != null, "Test \\Delta V_ext: Calculated v["+t+"] was null");
			
			// Compare to zero external influence case:
			RealVector v_expected_base = v_base[t].add(deltaV_ext[t]);
			for(int i = 0; i < 6; ++i)
			{
				assertTrue(v[t].getEntry(i) == v_expected_base.getEntry(i), "Test \\Delta V_ext: v["+t+"]_"+i+" not as expected; was "+v[t].getEntry(i)+", expected "+v_expected_base.getEntry(i));
			}
			
			// Compare with v_t case:
			RealVector v_expected = v_t.add(deltaV_ext[t]);
			for(int i = 0; i < 3; ++i)
			{
				assertTrue(equals(v[t].getEntry(i), v_expected.getEntry(i), 5e-4), "Test \\Delta V_ext: v["+t+"]_"+i+" not close to v_t+v_ext[t]; was "+v[t].getEntry(i)+", expected "+v_expected.getEntry(i));
			}
			for(int i = 3; i < 6; ++i)
			{
				assertTrue(equals(v[t].getEntry(i), v_expected.getEntry(i), 5e-3), "Test \\Delta V_ext: v["+t+"]_"+i+" not close to v_t+v_ext[t]; was "+v[t].getEntry(i)+", expected "+v_expected.getEntry(i));
			}
		}

		endTest("test_external()");
	}
}