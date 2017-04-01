package ellipsis.energy.casestudy;

import static ellipsis.energy.test.GridlessADPIEEE13BusGridTest.makeDG;
import static ellipsis.energy.test.IEEE13BusGrid.BASE_POWER;
import static ellipsis.energy.test.IEEE13BusGrid.BASE_VOLTAGE;

import javax.swing.JFrame;
import javax.swing.JOptionPane;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.DistributedSource;
import ellipsis.energy.grid.GridDisplay;
import ellipsis.energy.grid.Transformer;
import ellipsis.energy.test.IEEE13BusGrid;

public class CS05B_SensitivityAnalysis
{
	public static void main(String[] args)
	{
		IEEE13BusGrid grid = new IEEE13BusGrid();
		GridDisplay gridDisplay = new GridDisplay(grid, BASE_POWER, BASE_VOLTAGE);
		JFrame frame = gridDisplay.showInFrame();
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		grid.bus680.addChild(makeDG("DG680"));
		grid.bus675.addChild(makeDG("DG675"));		
		grid.bus684.addChild(makeDG("DG684"));
		grid.bus645.addChild(makeDG("DG645"));
		grid.bus611.addChild(makeDG("DG611"));
		grid.bus646.addChild(makeDG("DG646"));
		
//		testDifferentTapPositions(grid);
		
		testShiftInPower(grid, gridDisplay);
	}

	public static void testShiftInPower(IEEE13BusGrid grid, GridDisplay gridDisplay)
	{
		Transformer transformer = grid.getTransformer("650-");
		transformer.setT(1.0);
		
		// Analyse at minimum power output:
		System.out.println("Min DG output:");
		setDGMin(grid, "DG680");
		setDGMin(grid, "DG675");
		setDGMin(grid, "DG684");
		setDGMin(grid, "DG645");
		setDGMin(grid, "DG611");
		setDGMin(grid, "DG646");
		AnalysisResults results = analyse(grid);
		showMatrix(results.sensitivities());
		gridDisplay.analyse();
		JOptionPane.showMessageDialog(null, "Min DG output set.");
		
		// Analyse at maximum power output:
		System.out.println("Max DG output:");
		setDGMax(grid, "DG680");
		setDGMax(grid, "DG675");
		setDGMax(grid, "DG684");
		setDGMax(grid, "DG645");
		setDGMax(grid, "DG611");
		setDGMax(grid, "DG646");
		results = analyse(grid);
		showMatrix(results.sensitivities());
		gridDisplay.analyse();
		JOptionPane.showMessageDialog(null, "Max DG output set.");
	}

	public static void setDGMin(IEEE13BusGrid grid, String name)
	{
		DistributedSource dg = grid.getDistributedSource(name);
		dg.setPowerOutput(dg.getPmin(), 0);
	}

	public static void setDGMax(IEEE13BusGrid grid, String name)
	{
		DistributedSource dg = grid.getDistributedSource(name);
		dg.setPowerOutput(dg.getPmax(), 0);
	}

	public static void testDifferentTapPositions(IEEE13BusGrid grid)
	{
		Transformer transformer = grid.getTransformer("650-");
		RealMatrix sensitivities[] = new RealMatrix[11];
		testTapPositions(grid, transformer, sensitivities);
		logSensitivities(sensitivities);
	}

	private static void logSensitivities(RealMatrix[] sensitivities)
	{
		RealMatrix min = new Array2DRowRealMatrix(sensitivities[0].getRowDimension(), sensitivities[0].getColumnDimension());
		RealMatrix max = new Array2DRowRealMatrix(sensitivities[0].getRowDimension(), sensitivities[0].getColumnDimension());
		for(int i = 0; i < sensitivities.length; ++i)
		{
			System.out.println("m="+m(i));
			for(int r = 0; r < sensitivities[i].getRowDimension(); ++r)
			{
				for(int c = 0; c < sensitivities[i].getColumnDimension(); ++c)
				{
					double s_rc = sensitivities[i].getEntry(r, c);
					System.out.print(s_rc);
					System.out.print(',');
					
					double min_rc = min.getEntry(r, c);
					double max_rc = max.getEntry(r, c);
					if(min_rc == 0 || s_rc < min_rc)
						min.setEntry(r, c, s_rc);
					if(max_rc == 0 || s_rc > max_rc)
						max.setEntry(r, c, s_rc);
				}
				System.out.println();
			}
			System.out.println();
		}
		
		// Min:
		System.out.println("Min:");
		showMatrix(min);
		System.out.println("Max:");
		showMatrix(max);
		System.out.println("Range:");
		showMatrix(max.subtract(min));
	}

	public static void showMatrix(RealMatrix m)
	{
		for(int r = 0; r < m.getRowDimension(); ++r)
		{
			for(int c = 0; c < m.getColumnDimension(); ++c)
			{
				double m_rc = m.getEntry(r, c);
				System.out.print(m_rc);
				System.out.print(',');
			}
			System.out.println();
		}
	}

	public static void testTapPositions(IEEE13BusGrid grid, Transformer transformer, RealMatrix[] sensitivities)
	{
		for(int i = 0; i < 11; ++i)
		{
			double m = m(i);
			transformer.setT(m);
			AnalysisResults results = analyse(grid);
			sensitivities[i] = results.sensitivities();
		}
	}

	public static double m(int i)
	{
		return 0.95 + i*0.01;
	}

	public static AnalysisResults analyse(IEEE13BusGrid grid)
	{
		LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(BASE_POWER);
        lfa.setBaseVoltage(BASE_VOLTAGE);
        AnalysisResults results = lfa.analyse();
        if(!results.getDidConverge())
        	throw new RuntimeException("Alanysis did not converge!");
        return results;
	}
}