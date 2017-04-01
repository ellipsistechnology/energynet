package ellipsis.energy.test;

import javax.swing.JFrame;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.grid.Capacitor;
import ellipsis.energy.grid.DistributedSource;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.GridDisplay;
import ellipsis.energy.grid.Load;
import ellipsis.energy.grid.Source;
import ellipsis.energy.smartgrid.ControllableDemand;
import ellipsis.energy.smartgrid.ControllableDemand.CDLoad;

/**
 * Ref IEEE 06912981
 * @author bmillar
 *
 */
public class ConsensusTest
{
	public static void main(String[] args)
	{
		new ConsensusTest().run();
	}
	
	private Grid grid;

	private void run()
	{
		this.grid = makeGrid();
//		int i = 0;
//		System.out.println("Storage:");
//		for (ControllableDemand s : grid.get(ControllableDemand.class))
//		{
//			System.out.println(i+":"+s.getName());
//			++i;
//		}
		init();
		loop();
		update();
//		analyse();
		System.out.println("Final mismatch is "+mismatch());
	}

	private void update()
	{
		int i = 0;
		for (ControllableDemand s : grid.get(ControllableDemand.class))
		{
			s.setChargeRate(storagePower[i], 0);
			++i;
		}
	}

//	private void analyse()
//	{
//		GridDisplay.showInFrame(grid, IEEE13BusGrid.BASE_POWER, IEEE13BusGrid.BASE_VOLTAGE).setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//	}
	
	
	
	////  ////

	private static final int storageCount = 6;
	private static final double epsilon = 0.005;
	private static final int K = 100;
	private static final double[] loss_a = new double[]{0.85, 0.80, 0.83, 0.79, 0.92, 0.87};
	private static final double[] loss_b = new double[]{0.0077, 0.0082, 0.0065, 0.0075, 0.0085, 0.007};
	private static final double[][] d = new double[][]{
		{0.314,	0.4,	0.285,	0,  	0,  	0},
		{0.4, 	0.314,	0.285,	0,  	0,  	0},
		{0.285,	0.285,	-0.19,	0.285,	0.333,	0},
		{0,  	0,  	0.285,	0.214,	0,  	0.5},
		{0,  	0,  	0.333,	0,  	0.666,	0},
		{0,  	0,  	0,  	0.5,	0,  	0.5},
	};
	
	double[] r = new double[storageCount];
	double[] mismatch = new double[storageCount];
	double[] storagePower = new double[storageCount];
	
	private void init()
	{
		int i = 0;
		for (ControllableDemand s : grid.get(ControllableDemand.class))
		{
			storagePower[i] = s.getChargeRate().abs();
			mismatch[i] = mismatch();
			r[i] = loss_a[i] - 2*loss_b[i]*storagePower[i];
			++i;
		}
		log(0);
	}
	
	public void loop()
	{
		for(int k = 0; k < K; ++k)
		{
			step();
			log(k+1);
		}
	}
	
	private void log(int k)
	{
		System.out.print(k);
		System.out.print(',');
		for(int i = 0; i < storageCount; ++i)
		{
			System.out.print(r[i]);
			System.out.print(',');
		}
		for(int i = 0; i < storageCount; ++i)
		{
			System.out.print(mismatch[i]);
			System.out.print(',');
		}
		for(int i = 0; i < storageCount; ++i)
		{
			System.out.print(storagePower[i]);
			System.out.print(',');
		}
		System.out.println();
	}

	private void step()
	{
		// New values:
		double[] r_new = new double[storageCount];
		double[] storagePower_new = new double[storageCount];
		double[] mismatch_new = new double[storageCount];
		
		for(int i = 0; i < storageCount; ++i)
		{
			r_new[i] = nextIncrementalCost(i);
			storagePower_new[i] = nextStoragePower(i, r_new[i]);
		}
		
		for(int i = 0; i < storageCount; ++i)
		{
			mismatch_new[i] = nextMismatch(i, storagePower_new);
		}
		
		// Update current values:
		r = r_new;
		mismatch = mismatch_new;
		storagePower = storagePower_new;
	}

	private double nextMismatch(int i, double[] storagePower_new)
	{
		double mismatch_new = 0;
		for(int j = 0; j < storageCount; ++j)
		{
			mismatch_new += d[i][j]*(mismatch[j]+(storagePower_new[j]-storagePower[j]));
		}
		return mismatch_new;
	}

	private double nextStoragePower(int i, double r_i)
	{
		return -(r_i+loss_a[i])/(2*loss_b[i]);
	}

	private double nextIncrementalCost(int i)
	{
		double r_new = epsilon*mismatch[i];
		for(int j = 0; j < storageCount; ++j)
		{
			r_new += d[i][j]*r[j];
		}
		return r_new;
	}

	public double mismatch()
	{
		return totalGeneration() - totalLoad();
	}

	public double totalGeneration()
	{
		double p = 0;
		for (Source s : grid.getSources())
		{
			if(s instanceof DistributedSource)
			{
				p += s.getPowerOutput().abs();
			}
		}
		return p;
	}
	
	public double totalLoad()
	{
		double p = 0;
		for (Load l : grid.getLoads())
		{
			if(!(l instanceof CDLoad) && !(l instanceof Capacitor))
			{
				p += l.getLoad().abs();
			}
		}
		return p;
	}
	
	
	//// Setup ////

	private Grid makeGrid()
	{
		IEEE13BusGrid grid = new IEEE13BusGrid();
		
		grid.bus680.addChild(makeStorage("S680"));
		grid.bus675.addChild(makeStorage("S675"));		
		grid.bus684.addChild(makeStorage("S684"));
		grid.bus645.addChild(makeStorage("S645"));
		grid.bus611.addChild(makeStorage("S611"));
		grid.bus646.addChild(makeStorage("S646"));
		
		return grid;
	}

	public static DistributedSource makeDG(String name)
	{
		DistributedSource dg = new DistributedSource();
		dg.setName(name);
		dg.setPmax(1400e3);
		dg.setPowerOutput(1400e3, 0.0);
		return dg;
	}

	public static ControllableDemand makeStorage(String name)
	{
		ControllableDemand s1 = new ControllableDemand();
		s1.setName(name);
		s1.setMaxCapacity(10000e3);
		s1.setCapacity(5000e3);
		s1.setMaxChargeRate(5000e3);
		s1.setMaxDischargeRate(5000e3);
		s1.setChargeRate(new Complex(0.0));
		return s1;
	}
}