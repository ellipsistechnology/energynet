package ellipsis.energy.casestudy;

import static ellipsis.energy.test.GridlessADPIEEE13BusGridTest.makeGridlessADP;
import static ellipsis.energy.test.GridlessADPIEEE13BusGridTest.makeStorage;

import java.util.List;
import java.util.Random;

import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.test.GridTestHelper;
import ellipsis.energy.test.IEEE13BusGrid;
import ellipsis.energy.util.DiscreteGridlessADP;
import ellipsis.energy.util.GridlessADP;
import ellipsis.energy.util.GridlessADPFactory;
import ellipsis.util.Pair;

/**
 * Uses CS01's cost function.
 * Test cases compare (all full IEEE 13 bus network):
 * 	DP,
 * 	Random,
 *  ADP (non-exploitative),
 *  ADP (exploitative),
 *  ADP (non-exploitative and successive).
 *  
 * @author bmillar
 *
 */
public class CS02 extends GridTestHelper
{
	public static void main(String[] args)
	{
		assertAssertsOn();
		CS02 cs = new CS02();
		cs.parseArgs(args);
		cs.run(0);
		
		// Average results:
		if(cs.iterations > 0)
		{
			for(int i = 0; i < cs.iterations; ++i)
			{
				cs.run(i);
			}
			for(int t = 0; t < cs.plainAverageCtgs.length; ++t)
			{
				cs.plainAverageCtgs[t] /= cs.iterations+1;
				cs.stadpAverageCtgs[t] /= cs.iterations+1;
				cs.estadpAverageCtgs[t] /= cs.iterations+1;
			}
			
			System.out.println("Average:\nt,plainCtg,stadpCtg,estadpCtg");
			for(int t = 0; t < cs.plainAverageCtgs.length; ++t)
			{
				double plainCtg = cs.plainAverageCtgs[t];
				double stadpCtg = cs.stadpAverageCtgs[t];
				double estadpCtg = cs.estadpAverageCtgs[t];
				System.out.println(t+","+plainCtg+","+stadpCtg+","+estadpCtg);
			}
		}
	}

	private double[] plainAverageCtgs = new double[13];
	private double[] stadpAverageCtgs = new double[13];
	private double[] estadpAverageCtgs =  new double[13];
	public void run(int seed)
	{
		// Setup ADP:
		GridlessADP estadp = makeGrid();
		estadp.setNr(20);
		estadp.setCostFunction(new CS01.CS01CostFunction());
//				new CostFunction() 
//		{
//			@Override
//			public double g(GridlessADP adp, int t, RealVector x, RealVector u)
//			{
//				RealVector v = adp.voltages(t, x);
////				RealVector PQ = adp.powerFromState(x);
//				
//				int dimension = v.getDimension();
//				
//				RealVector vCostVector = v.subtract(new ArrayRealVector(dimension, 1.0));
//				double vCost = vCostVector.dotProduct(vCostVector);
//				
//				return vCost;
//			}
//		};
		
		// Deterministic comparison:
		if(verbose)
			System.out.println("DP...");
		String dpScheduleFileName = "/opt/energynet/cs02/dpSchedule";
		@SuppressWarnings({ "unchecked", "rawtypes" })
		List<Pair<RealVector, Double>> dpSchedule = 
			(List<Pair<RealVector, Double>>) convertFromSerializable((List<Pair>) read(dpScheduleFileName));
		if(dpSchedule == null)
			dpSchedule = dp((DiscreteGridlessADP)estadp);
		save(convertToSerializable(dpSchedule), dpScheduleFileName);
		
		// Random:
		if(verbose)
			System.out.println("Random...");
		String randomCtgsFileName = "/opt/energynet/cs02/randomSchedule";
		double[] randomCtgs = (double[]) read(randomCtgsFileName);
		if(randomCtgs == null)
		{
			randomCtgs = new double[estadp.getT()+1];
			int randomIterations = 1000;
			Random rand = new Random(0);
			for(int i = 0; i < randomIterations; ++i)
			{
				double[] ctgs = new double[estadp.getT()+1];
				randomPath((DiscreteGridlessADP)estadp, ctgs, rand);
				for (int j = 0; j < ctgs.length; j++)
				{
					randomCtgs[j] += ctgs[j];
				}
			}
			for (int j = 0; j < randomCtgs.length; j++)
			{
				randomCtgs[j] /= randomIterations;
			}
			save(randomCtgs, randomCtgsFileName);
		}
		
		// ADP:
		if(verbose)
			System.out.println("ADP...");
		double[][] h = new double[][]{{0.5}};//[estadp.storageCount];//+2*(estadp.dgCount+estadp.storageCount+estadp.loadCount)];
//		h[0] = 0.5; // SOC
//		for(int i = 1; i < h.length; ++i)
//			h[i] = 0.5; // DG, storage, loads
		
		DiscreteGridlessADP adp = GridlessADPFactory.copyInto((DiscreteGridlessADP)estadp, new CS02ADP(false, true));
		adp.rand = new Random(seed);
		adp.train(h);
		
		double[] adpCtgs = new double[estadp.getT()+1];
		/*RealVector[] stadpSchedule = */adp.schedule(adpCtgs);
//		List<Pair<RealVector, Double>> adpCtgs_0 = ctg_0(adp);
		
		// EST-ADP (Exploitative):
		if(verbose)
			System.out.println("EST-ADP...");

		((DiscreteGridlessADP)estadp).rand = new Random(seed);
		estadp.train(h);
		
		double[] estadpCtgs = new double[estadp.getT()+1];
		/*RealVector[] estadpSchedule = */estadp.schedule(estadpCtgs);
//		List<Pair<RealVector, Double>> estadpCtgs_0 = ctg_0(estadp);
		
		// ST-ADP:
		if(verbose)
			System.out.println("ST-ADP...");
		DiscreteGridlessADP stadp = GridlessADPFactory.copyInto((DiscreteGridlessADP)estadp, new CS02ADP(true, false));

		stadp.rand = new Random(seed);
		stadp.train(h);
		
		double[] stadpCtgs = new double[estadp.getT()+1];
		/*RealVector[] stadpSchedule = */stadp.schedule(stadpCtgs);
//		List<Pair<RealVector, Double>> stadpCtgs_0 = ctg_0(stadp);
		
		// Print comparison:
		System.out.println("t,randomCtg,dpCtg,plainCtg,stadpCtg,estadpCtg");
		for(int t = 0; t <= estadp.getT(); ++t)
		{
			double randomCtg = randomCtgs[t];
			double dpCtg = dpSchedule.get(t).getValue();
			double plainCtg = adpCtgs[t];
			double stadpCtg = stadpCtgs[t];
			double estadpCtg = estadpCtgs[t];
			System.out.println(t+","+randomCtg+","+dpCtg+","+plainCtg+","+stadpCtg+","+estadpCtg);
			
			plainAverageCtgs[t] += plainCtg;
			stadpAverageCtgs[t] += stadpCtg;
			estadpAverageCtgs[t] += estadpCtg;
		}
		
//		System.out.println("u,plain,st,est");
//		for (int i = 0; i < adpCtgs_0.size(); i++)
//		{
//			Pair<RealVector, Double> plain = adpCtgs_0.get(i);
//			RealVector plain_u_0 = plain.getKey();
//			double plain_ctg_0 = plain.getValue();
//			
//			Pair<RealVector, Double> st = stadpCtgs_0.get(i);
//			RealVector st_u_0 = st.getKey();
//			double st_ctg_0 = st.getValue();
//			
//			Pair<RealVector, Double> est = estadpCtgs_0.get(i);
//			RealVector est_u_0 = est.getKey();
//			double est_ctg_0 = est.getValue();
//			
//			assertTrue(plain_u_0.equals(est_u_0) && plain_u_0.equals(st_u_0), "Controls do not match: "+plain_u_0+","+st_u_0+","+est_u_0);
//			
//			System.out.println(VectorHelper.printVector(plain_u_0)+plain_ctg_0+","+st_ctg_0+","+est_ctg_0);
//		}
	}

//	private List<Pair<RealVector, Double>> ctg_0(GridlessADP adp)
//	{
//		List<Pair<RealVector, Double>> points = new ArrayList<>();
//		Set<RealVector> U = adp.enumerateControls(0, adp.x_0);
//		for (RealVector u : U)
//		{
//			RealVector x_1 = adp.f(adp.x_0, u, adp.wZero(0));
//			double ctg_1 = adp.approxJ[0].value(x_1);
//			points.add(new Pair<RealVector, Double>(u, ctg_1));
//		}
//		return points;
//	}

	public GridlessADP makeGrid()
	{
		IEEE13BusGrid grid = new IEEE13BusGrid();
		
		grid.bus632.addChild(makeStorage("S1"));
//		grid.bus671.addChild(makeStorage("S2"));
//
//		grid.bus680.addChild(makeDG("DG1"));
//		grid.bus675.addChild(makeDG("DG2"));		
//		grid.bus684.addChild(makeDG("DG3"));
//		grid.bus645.addChild(makeDG("DG4"));
//		grid.bus611.addChild(makeDG("DG5"));
//		grid.bus646.addChild(makeDG("DG6"));
//		grid.bus634.addChild(makeDG("DG7"));
//		grid.bus671.addChild(makeDG("DG8"));
		
		GridlessADP adp = makeGridlessADP(grid);
		
		return adp;
	}
	
	public static class CS02ADP extends DiscreteGridlessADP
	{
		private boolean successive;
		private boolean pureExploitative;

		public CS02ADP(boolean successive, boolean pureExploitative)
		{
			this.successive = successive;
			this.pureExploitative = pureExploitative;
		}

		public double explorationRate(int k, int t)
		{
			if(!pureExploitative)
				return 
					k < N_r*(T-t-1) ? 1 : // all random before training
					k > N_r*(T-t)   ? 0 : // only optimal after training
					                  1 ; // random during training
			else
					return 0;
		}
		
		@Override
		public boolean trained(int k, int t)
		{
			if(successive)
				return super.trained(k, t);
			else
				return false;
		}
		
		/*@Override
		public void setupEstimator(double[] defaultBandwidth, NonParametricApproximator estimator)
		{
			super.setupEstimator(defaultBandwidth, estimator);
			estimator.setCorrectInterpolation(false);
		}*/
	}
}