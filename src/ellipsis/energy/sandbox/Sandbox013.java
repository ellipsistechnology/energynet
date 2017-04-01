package ellipsis.energy.sandbox;

import static ellipsis.energy.sandbox.GridConstants.BASE_IMPEDANCE;
import static ellipsis.energy.sandbox.GridConstants.BASE_POWER;
import static ellipsis.energy.sandbox.GridConstants.BASE_VOLTAGE;
import static ellipsis.energy.sandbox.GridConstants.BUS_2;
import static ellipsis.energy.sandbox.GridConstants.BUS_3;
import static ellipsis.energy.sandbox.GridConstants.BUS_4;
import static ellipsis.energy.sandbox.GridConstants.BUS_5;
import static ellipsis.energy.sandbox.GridConstants.DG_3;
import static ellipsis.energy.sandbox.GridConstants.DG_4;
import static ellipsis.energy.sandbox.GridConstants.LINE_1_2;
import static ellipsis.energy.sandbox.GridConstants.LINE_1_3;
import static ellipsis.energy.sandbox.GridConstants.LINE_2_3;
import static ellipsis.energy.sandbox.GridConstants.LINE_2_4;
import static ellipsis.energy.sandbox.GridConstants.LINE_3_5;
import static ellipsis.energy.sandbox.GridConstants.LOAD_2;
import static ellipsis.energy.sandbox.GridConstants.LOAD_5;
import static ellipsis.energy.sandbox.GridConstants.SLACK_BUS;
import static ellipsis.energy.sandbox.GridConstants.SLACK_SOURCE;
import static ellipsis.energy.sandbox.Sandbox004.maxPower;
import static ellipsis.energy.sandbox.Sandbox004.p0;
import static ellipsis.energy.sandbox.Sandbox004.q0;
import static ellipsis.energy.sandbox.Sandbox008.debugLogNames;
import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Source;

/**
 * Min/max Lagrange gradient method for OPF.
 * Uses real and imaginary parts of voltage.
 * Uses augmented Lagrangian.
 * @author bmillar
 *
 */
public class Sandbox013
{	
	public static void main(String[] args)
	{
		new Sandbox013().run();
	}
	
	public static interface IndexedFunction
	{
		double value(int i);
	}
	
	public static interface IndexedBooleanFunction
	{
		boolean value(int i);
	}
	
//	private static final double errorTolerance = 1e-3; // Maximum value of g(x) that's considered close enough to zero.
//	private static final double minConstraintImprovement = 0.9; // Minimum improvement percentage of g(x) constraint
//	private static final double MAX_AUG_SCALE = 1e12;
//	private static final double AUG_SCALE_STEP = 2;
	
	private static final double V_MAX = 1.05;
	private static final double V_MIN = 0.95;
	private static final double CON_ERROR_TOLERANCE = 1e-3; // Maximum value of g(x) that's considered close enough to zero.
	private static final double MIN_CONSTRAINT_IMPROVEMENT = 0.9; // Minimum improvement percentage of g(x) constraint
	private static final double MAX_AUG_SCALE = 1e0;
	private static final double AUG_SCALE_STEP = 1.1;
	private static final int K = 1000;
	
	private static boolean ZERO_SMALL_G = false;
	private static boolean PROJECT_X = true;
	
	protected double gamma = 0.1;
	protected double c = 1e-3;//0.001;
	protected RealVector p_max;
	protected RealVector generatorMask;
	
	protected RealVector p, q, e, f, lp, lq;
	protected RealMatrix G, B;
	private int slackIndex;
	private Grid grid;
	
	protected void run()
	{
		AnalysisResults results = init();
		loop();
		testResults(results);
	}
	
	protected void testResults(AnalysisResults results)
	{
		// Set values to results:
		Map<String, Integer> busNumbers = results.getBusNumbers();
		for (String busName : busNumbers.keySet())
		{
			Bus bus = grid.getBus(busName);
			Source dg = null;
			for (Source source : bus.getSources())
			{
				dg = source;
				break;
			}
			if(dg != null)
			{
				int i = busNumbers.get(busName);
				dg.setPowerOutput(BASE_POWER*p(i), BASE_POWER*q(i));
			}
		}
//        grid.getSource(DG_3).setPowerOutput(BASE_POWER*6.5, BASE_POWER*0);
//        grid.getSource(DG_4).setPowerOutput(BASE_POWER*6.5, BASE_POWER*0);
		
		// Analyse grid to get parameters:
		LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(grid.getBasePower());
        lfa.setBaseVoltage(grid.getBaseVoltage());
        AnalysisResults results2 = lfa.analyse(complexVoltages());
        
        System.out.println("\nFinal Values:");
        System.out.println("index,bus,e,f,p,q,,|v|,|s|");
        for (String busName : busNumbers.keySet())
		{
			int i = busNumbers.get(busName);
        	double e_i = e(i);
			double f_i = f(i);
			String v_abs = format(Math.sqrt(e_i*e_i+f_i*f_i));
			double p_i = p(i);
			double q_i = q(i);
			String s_abs = format(Math.sqrt(p_i*p_i+q_i*q_i));
			System.out.println(i+","+busName+','+format(e_i)+','+format(f_i)+','+format(p_i)+','+format(q_i)+",,"+v_abs+','+s_abs);
        }
        System.out.println("\nAnalysis Results:");
        logBusses(results2);
	}

	private Complex[] complexVoltages()
	{
		int dimension = e.getDimension();
		Complex[] v = new Complex[dimension];
		for(int i = 0; i < dimension; ++i)
		{
			v[i] = new Complex(e(i), f(i));
		}
		return v;
	}

	protected AnalysisResults init()
	{
		initGrid();
		
		// Analyse grid to get parameters:
		LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(grid.getBasePower());
        lfa.setBaseVoltage(grid.getBaseVoltage());
        AnalysisResults results = lfa.analyse();
        
        logBusses(results);
        
        // Initialise variables:
        initVariables(results);
        
        // Debug:
		debugHeader(results);
		
		return results;
	}
	
	private void initGrid()
	{
		// Copied from Sandbox008:
		grid = 
            Grid.grid().
            
                Bus(SLACK_BUS).
                
                    SlackSource(SLACK_SOURCE, 1.00*BASE_VOLTAGE, 0, 0, 0).
                    
                    Line(LINE_1_2, 1, 0.02*BASE_IMPEDANCE, 0.04*BASE_IMPEDANCE).
                        Bus(BUS_2).
                            Load(LOAD_2).
                            Line(LINE_2_3, BUS_3, 1, 0.0125*BASE_IMPEDANCE, 0.025*BASE_IMPEDANCE).
                            Line(LINE_2_4, 1.5, 0.0125*BASE_IMPEDANCE, 0.025*BASE_IMPEDANCE).
	                        	Bus(BUS_4).
	                        		DistributedSource(DG_4).
	                        	terminate().
	                          terminate().
                        terminate().
                    terminate().
                    
                    Line(LINE_1_3, 1, 0.01*BASE_IMPEDANCE, 0.03*BASE_IMPEDANCE).
                        Bus(BUS_3).
                        	DistributedSource(DG_3).
                        	Line(LINE_3_5, 1.5, 0.0125*BASE_IMPEDANCE, 0.025*BASE_IMPEDANCE).
	                        	Bus(BUS_5).
	                        		Load(LOAD_5).
	                        	terminate().
	                        terminate().
                        terminate().
                    terminate().
                    
                terminate().
                
            grid();
		
		grid.setBaseVoltage(BASE_VOLTAGE);
		grid.setBasePower(BASE_POWER);
		
		grid.getLoad(LOAD_2).setLoad(new Complex(150e6, 50e6));
		grid.getLoad(LOAD_5).setLoad(new Complex(80e6, 25e6));
		
        grid.getSource(DG_3).setPmax(650e6);
        grid.getSource(DG_3).setPowerOutput(100e6, 0);
        grid.getSource(DG_4).setPmax(650e6);
        grid.getSource(DG_4).setPowerOutput(100e6, 0);
	}

	protected void initVariables(AnalysisResults results)
	{
        slackIndex = results.getSlackIndex();
        
		// Get admittances:
        FieldMatrix<Complex> Y = results.getAdmittanceMatrix().Y;
        int dimension = Y.getRowDimension();
        B = new Array2DRowRealMatrix(dimension, dimension);
        G = new Array2DRowRealMatrix(dimension, dimension);
        for(int i = 0; i < dimension; ++i)
        {
        	for(int j = 0; j < dimension; ++j)
        	{
        		Complex Y_ij = Y.getEntry(i, j);
				G.setEntry(i, j, Y_ij.getReal());
				B.setEntry(i, j, Y_ij.getImaginary());
        	}
        }
        
        // Prepare generator mask (1 if generator, 0 otherwise):
        generatorMask = new ArrayRealVector(dimension);
        for (String name : results.getBusNumbers().keySet())
		{
			int index = results.getBusNumbers().get(name);
			Bus bus = grid.getBus(name);
			if(!bus.getGeneratedPower().equals(Complex.ZERO))
				generatorMask.setEntry(index, 1);
		}
        
        // Init power vector:
        p = p0(results);
        q = q0(results);
        p.setEntry(slackIndex, 0);
        q.setEntry(slackIndex, 0);
        
        p_max = maxPower(results, grid);
        
        // Init voltage vectors:
        e = new ArrayRealVector(dimension, 1);
        f = new ArrayRealVector(dimension, 0);
        
        // Init Lagrange multipliers:
        lp = new ArrayRealVector(dimension, 0);
        lq = new ArrayRealVector(dimension, 0);
	}
	
	public void loop()
	{
		int dimension = p.getDimension();
		IndexedBooleanFunction generators = i -> generatorMask.getEntry(i) == 1;
		IndexedBooleanFunction nonSlack = i -> i != slackIndex;

		for(int k = 0; k < K; ++k)
		{
			step(dimension, generators, nonSlack, k);
		}
	}

	private boolean debug_backtrack = false;
	private void step(int dimension, IndexedBooleanFunction generators, IndexedBooleanFunction nonSlack, int k)
	{
		// Remember current constraint values (before x step):
		RealVector gp = gp(p, q, e, f);
		RealVector gq = gq(p, q, e, f);

		// Step:
		double _t = stepX(dimension, generators, nonSlack);
		stepLambda(dimension, nonSlack);
		stepAugScale(gp, gq);

		// Log:
		debug(k, _t); // t has been updated at end of loop so un-update it.
	}

	protected double stepX(int dimension, IndexedBooleanFunction generators, IndexedBooleanFunction nonSlack)
	{
		// Find x step direction:
		RealVector dp = mask(vectorByComponents(this::gradL_p, dimension), generators);
		RealVector dq = mask(vectorByComponents(this::gradL_q, dimension), generators);
		RealVector de = mask(vectorByComponents(this::gradL_e, dimension), nonSlack);
		RealVector df = mask(vectorByComponents(this::gradL_f, dimension), nonSlack);
			
			// Check stopping criterion:
	//			double gradNorm = dp.getL1Norm() + dq.getL1Norm() + de.getL1Norm() + df.getL1Norm();
	//			if(gradNorm < 10e-6)
	//				break;
			
		// Backtrack to find best step size:
		double t = 1;//0.1; // Step size
		double alpha = 0.5; // Improvement scale
		double beta = 0.5; // Step update size.
		double L = lagrangian(p, q, e, f, lp, lq);
		if(debug_backtrack)
			checkGradient(dp, dq, de, df);
		double L_check = L();
		double L_new;
		RealVector p_new;
		RealVector q_new;
		RealVector e_new;
		RealVector f_new;
		double step = alpha*(dp.dotProduct(dp) + dq.dotProduct(dq) + de.dotProduct(de) + df.dotProduct(df));
		double minImprovement;
	
		do
		{
			p_new = p.subtract(dp.mapMultiply(t));
			q_new = q.subtract(dq.mapMultiply(t));
			e_new = e.subtract(de.mapMultiply(t));
			f_new = f.subtract(df.mapMultiply(t));
			
			L_new = lagrangian(p_new, q_new, e_new, f_new, lp, lq);
			minImprovement = L - t*step;
			
			if(debug_backtrack)
			{
				double lowerBound = L - t*step/alpha;
				System.out.println(t+","+L_new+","+minImprovement+","+lowerBound);
			}
			
			t = t*beta;
		} while(L_new > minImprovement);
		
		if(debug_backtrack)
		{
			double lowerBound = L; // t = 0
			System.out.println("0,"+L+","+L+","+lowerBound);
			L_new = lagrangian(
				p.subtract(dp.mapMultiply(-0.01)), 
				q.subtract(dq.mapMultiply(-0.01)), 
				e.subtract(de.mapMultiply(-0.01)),
				f.subtract(df.mapMultiply(-0.01)), 
				lp, lq);
			lowerBound = L - (-0.01)*step/alpha; // t = -0.01
			System.out.println("-0.01,"+L_new+","+(L-(-0.01)*step)+","+lowerBound);
		}
		
		// Project onto constrained set X:
		if(PROJECT_X)
		{
			projectVoltage(e_new, f_new);
			projectPowers(p_new, q_new); // FIXME only for generators
		}
		
		// Update x: FIXME don't update P and Q for generators
		p = p_new;
		q = q_new;
		e = e_new;
		f = f_new;
		
		if(debug_backtrack)
		{
			L_check = L();
			L_new = lagrangian(p, q, e, f, lp, lq);
			System.out.println(L_check+" ?= "+L_new);
		}
		
		double _t = t/beta;
		
		return _t;
	}

	protected void stepLambda(int dimension, IndexedBooleanFunction nonSlack)
	{
		// Maximise (l := l - c*g(x), ref http://en.wikipedia.org/wiki/Augmented_Lagrangian_method
		//           and Bertsekas' book):
		RealVector dlp = mask(vectorByComponents(this::gp/*this::gradL_lp*/, dimension), nonSlack);
		lp = lp.add(dlp.mapMultiply(c));
		RealVector dlq = mask(vectorByComponents(this::gq/*this::gradL_lq*/, dimension), nonSlack);
		lq = lq.add(dlq.mapMultiply(c));
		
		lp = project(lp, -1e14, 1e14);
		lq = project(lq, -1e14, 1e14);
	}

//	private boolean debug_augscale = false;
	protected void stepAugScale(RealVector gp, RealVector gq)
	{
		if(gp(p, q, e, f).getNorm() > MIN_CONSTRAINT_IMPROVEMENT*gp.getNorm() || gq(p, q, e, f).getNorm() > MIN_CONSTRAINT_IMPROVEMENT*gq.getNorm())
		{
			if(c < MAX_AUG_SCALE)
				c *= AUG_SCALE_STEP;
		}
		
//		if(debug_augscale)
//		{
//			System.out.println(gp(p, q, e, f).getNorm()+" ?> "+delta+"*"+gp.getNorm()+"="+delta*gp.getNorm());
//			System.out.println(gq(p, q, e, f).getNorm()+" ?> "+delta+"*"+gq.getNorm()+"="+delta*gq.getNorm());
//		}
	}

	private void projectVoltage(RealVector e, RealVector f)
	{
		int dimension = e.getDimension();
		for(int i = 0; i < dimension; ++i)
		{
			Complex v = new Complex(e.getEntry(i), f.getEntry(i));
			double abs = v.abs();
			if(abs > V_MAX)
			{
				v = v.multiply(V_MAX/abs);
			}
			else if(abs < V_MIN)
			{
				v = v.multiply(V_MIN/abs);
			}
			e.setEntry(i, v.getReal());
			f.setEntry(i, v.getImaginary());
		}
	}

	private void projectPowers(RealVector p, RealVector q)
	{
		int dimension = p.getDimension();
		for(int i = 0; i < dimension; ++i)
		{
			// Skip non-generators:
			if(generatorMask.getEntry(i) == 0)
				continue;
			
			// Restrict power magnitude:
			Complex s = new Complex(p.getEntry(i), q.getEntry(i));
			double abs = s.abs();
			double p_max_i = p_max.getEntry(i);
			if(abs > p_max_i)
			{
				s = s.multiply(p_max_i/abs);
			}
			p.setEntry(i, s.getReal());
			q.setEntry(i, s.getImaginary());
		}
	}

	private void checkGradient(RealVector dp, RealVector dq, RealVector de, RealVector df)
	{
		int dimension = p.getDimension();
		double delta = 1e-6;
		for(int i = 0; i < dimension; ++i)
		{
			// Estimate p slope:
			RealVector leftP = new ArrayRealVector(p);
			leftP.setEntry(i, p.getEntry(i)-delta);
			double left = lagrangian(leftP, q, e, f, lp, lq);
			
			RealVector rightP = new ArrayRealVector(p);
			rightP.setEntry(i, p.getEntry(i)+delta);
			double right = lagrangian(rightP, q, e, f, lp, lq);
			
			double pSlope = (right-left)/(2*delta);

			// Estimate p slope:
			RealVector leftQ = new ArrayRealVector(q);
			leftQ.setEntry(i, q.getEntry(i)-delta);
			left = lagrangian(p, leftQ, e, f, lp, lq);
			
			RealVector rightQ = new ArrayRealVector(q);
			rightP.setEntry(i, q.getEntry(i)+delta);
			right = lagrangian(p, rightQ, e, f, lp, lq);
			
			double qSlope = (right-left)/(2*delta);
			
			// Estimate p slope:
			RealVector leftE = new ArrayRealVector(e);
			leftE.setEntry(i, e.getEntry(i)-delta);
			left = lagrangian(p, q, leftE, f, lp, lq);
			
			RealVector rightE = new ArrayRealVector(e);
			rightE.setEntry(i, e.getEntry(i)+delta);
			right = lagrangian(p, q, rightE, f, lp, lq);
			
			double eSlope = (right-left)/(2*delta);

			// Estimate p slope:
			RealVector leftF = new ArrayRealVector(f);
			leftF.setEntry(i, f.getEntry(i)-delta);
			left = lagrangian(p, q, e, leftF, lp, lq);
			
			RealVector rightF = new ArrayRealVector(f);
			rightF.setEntry(i, f.getEntry(i)+delta);
			right = lagrangian(p, q, e, rightF, lp, lq);
			
			double fSlope = (right-left)/(2*delta);
			
			// Log:
			System.out.println(i
					+","+dp.getEntry(i)+","+pSlope
					+","+dq.getEntry(i)+","+qSlope
					+","+de.getEntry(i)+","+eSlope
					+","+df.getEntry(i)+","+fSlope);
		}
	}

	protected RealVector project(RealVector v, double low, double high)
	{
		int dimension = v.getDimension();
		RealVector vNew = new ArrayRealVector(dimension);
		for(int i = 0; i < dimension; ++i)
		{
			double v_i = v.getEntry(i);
			vNew.setEntry(i, project(v_i, low, high));
		}
		return vNew;
	}

	protected double project(double v, double low, double high)
	{
		return min(max(low, v), high);
	}

	private double lagrangian(RealVector p, RealVector q, RealVector e, RealVector f, RealVector lp, RealVector lq)
	{
		RealVector gp = gp(p, q, e, f);
		RealVector gq = gq(p, q, e, f);
		gp.setEntry(slackIndex, 0);
		gq.setEntry(slackIndex, 0);
		
		double aug = gp.dotProduct(gp) + gq.dotProduct(gq);
		double lgp = lp.dotProduct(gp);
		double lgq = lq.dotProduct(gq);
		
		return cost(p, q) + lgp + lgq + 0.5*c*aug;
	}

	protected double gradL_p(int i)
	{
		return 	gradC_p(i) - lp(i)

				// Penalty:
				-c*gp(i)
				;
	}
	
	/**
	 * Assuming c_i(p_i) = 0.5*(p_i - p_max)^2
	 * => dF_i/dp_i = p_i - p_max
	 * @param i
	 * @return
	 */
	private double gradC_p(int i)
	{
		return p(i) - p_max.getEntry(i);
	}

	protected double gradL_q(int i)
	{
//System.out.println(gradF_q(i));
//System.out.println(lq(i));
//System.out.println(c*gq(i));
		return gradC_q(i) - lq(i)

				// Penalty:
				-c*gq(i)
				;
	}
	
	/**
	 * Assuming c_i(q_i) = 0.5*q_i^2
	 * => dF_i/dp_i = q_i
	 * @param i
	 * @return
	 */
	private double gradC_q(int i)
	{
		return q(i);
	}

	protected double gradL_e(int j)
	{
		int dimension = p.getDimension();
		double constraintTerms = 
				 sum(i -> lp(i)*(e(i)*G(i,j) + f(i)*B(i,j)), dimension, j) // i != j
				+sum(i -> lq(i)*(f(i)*G(i,j) - e(i)*B(i,j)), dimension, j) // i != j
				+lp(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) // n != j
				        +2*G(j,j)*e(j))
				+lq(j)*(sum(n -> -f(n)*G(j,n) - e(n)*B(j,n), dimension, j) // n != j
				        -2*B(j,j)*e(j));
		
		double pt1 = c*sum(i -> i == slackIndex ? 0 : 
							 gp(i)*(e(i)*G(i,j) + f(i)*B(i,j))
							+gq(i)*(f(i)*G(i,j) - e(i)*B(i,j)),
							dimension, j);
		double pt2 = c*gp(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) 
		          +2*e(j)*G(j,j));
		double pt3 = c*gq(j)*(sum(n -> -f(n)*G(j,n) - e(n)*B(j,n), dimension, j) 
		          -2*e(j)*B(j,j));
		double penaltyTerms = pt1
							+
							(j == slackIndex ? 0 : 
							(
								+pt2
								+pt3
							));
		return
			constraintTerms
			+penaltyTerms
			;
	}
	
	private RealVector gp(RealVector p, RealVector q, RealVector e, RealVector f)
	{
		int dimension = p.getDimension();
		return vectorByComponents(i -> gp(p, q, e, f, i), dimension);
	}
	
	private double gp(RealVector p, RealVector q, RealVector e, RealVector f, int i)
	{
		if(i == slackIndex)
			return 0;
		
		int dimension = p.getDimension();
		double p_i = p.getEntry(i);
		double e_i = e.getEntry(i);
		double f_i = f.getEntry(i);
		double d = sum(n -> e_i*e.getEntry(n)*G(i,n) + f_i*f.getEntry(n)*G(i,n) + f_i*e.getEntry(n)*B(i,n) - e_i*f.getEntry(n)*B(i,n), dimension)
					- p_i;
		return ZERO_SMALL_G && (Math.abs(d) < CON_ERROR_TOLERANCE) ? 0 : d;
	}
	
	private RealVector gq(RealVector p, RealVector q, RealVector e, RealVector f)
	{
		int dimension = q.getDimension();
		return vectorByComponents(i -> gq(p, q, e, f, i), dimension);
	}
	
	private double gq(RealVector p, RealVector q, RealVector e, RealVector f, int i)
	{
		if(i == slackIndex)
			return 0;
		
		int dimension = p.getDimension();
		double q_i = q.getEntry(i);
		double e_i = e.getEntry(i);
		double f_i = f.getEntry(i);
		double d = sum(j -> f_i*e.getEntry(j)*G(i,j) - e_i*f.getEntry(j)*G(i,j) - e_i*e.getEntry(j)*B(i,j) - f_i*f.getEntry(j)*B(i,j), dimension)
					- q_i;
		return ZERO_SMALL_G && (Math.abs(d) < CON_ERROR_TOLERANCE) ? 0 : d;
	}

	private boolean debug_gp = false;
	private double gp(int i)
	{
		int dimension = p.getDimension();
		double sum = 0;
		for(int n = 0; n < dimension; ++n)
		{
			sum += e(i)*e(n)*G(i,n) + f(i)*f(n)*G(i,n) + f(i)*e(n)*B(i,n) - e(i)*f(n)*B(i,n);
			if(debug_gp)
			{
				System.out.println(e(i)+","+e(n)+","+G(i,n)+","+f(i)+","+f(n)+","+B(i,n));
			}
		}
		
		double d = sum//sum(n -> e(i)*e(n)*G(i,n) + f(i)*f(n)*G(i,n) + f(i)*e(n)*B(i,n) - e(i)*f(n)*B(i,n), dimension)
					- p(i);
		return ZERO_SMALL_G && (Math.abs(d) < CON_ERROR_TOLERANCE) ? 0 : d;
	}
	
	private double gq(int i)
	{
		int dimension = q.getDimension();
		double sum = sum(n -> f(i)*e(n)*G(i,n) - e(i)*f(n)*G(i,n) - e(i)*e(n)*B(i,n) - f(i)*f(n)*B(i,n), dimension);
		double d = sum - q(i);
		return ZERO_SMALL_G && (Math.abs(d) < CON_ERROR_TOLERANCE) ? 0 : d;
	}

	protected double gradL_f(int j)
	{
		int dimension = p.getDimension();
		double constraintTerms = sum(i -> lp(i)*(f(i)*G(i,j) - e(i)*B(i,j)), dimension, j) // i != j
					+sum(i -> lq(i)*(-e(i)*G(i,j) - f(i)*B(i,j)), dimension, j) // i != j
					+lp(j)*(sum(n -> f(n)*G(j,n) + e(n)*B(j,n), dimension, j) // n != j
					        +2*G(j,j)*f(j))
					+lq(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) // n != j
					        -2*B(j,j)*f(j));
		
		double pt1 = c*sum(i -> i == slackIndex ? 0 : gp(i)*( f(i)*G(i,j) - e(i)*B(i,j))
 			                                          +gq(i)*(-e(i)*G(i,j) - f(i)*B(i,j)),
 			   		       dimension, j);
		double pt2 = c*gp(j)*(sum(n -> f(n)*G(j,n) + e(n)*B(j,n), dimension, j) 
				              +2*f(j)*G(j,j));
		double pt3 = c*gq(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) 
				              -2*f(j)*B(j,j));
		double penaltyTerms = pt1
							+ 
							(j == slackIndex ? 0 : 
							(
								+pt2
								+pt3
							));
		
		return 
			constraintTerms
			+penaltyTerms
			;
	}
	
//	protected double gradL_lp(int i)
//	{
//		return 
//			sum(j -> e(i)*e(j)*G(i,j) + f(i)*f(j)*G(i,j) + f(i)*e(j)*B(i, j) - e(i)*f(j)*B(i,j), p.getDimension())
//			- p(i);
//	}
//	
//	protected double gradL_lq(int i)
//	{
//		return 
//			sum(j -> f(i)*e(j)*G(i,j) - e(i)*f(j)*G(i,j) - e(i)*e(j)*B(i, j) - f(i)*f(j)*B(i,j), p.getDimension())
//			- q(i);
//	}

	protected double p(int i)
	{
		return p.getEntry(i);
	}

	protected double q(int i)
	{
		return q.getEntry(i);
	}

	protected double lp(int i)
	{
		return lp.getEntry(i);
	}

	protected double lq(int i)
	{
		return lq.getEntry(i);
	}

	protected double B(int i, int j)
	{
		return B.getEntry(i, j);
	}

	protected double G(int i, int j)
	{
		return G.getEntry(i, j);
	}

	protected double e(int i)
	{
		return e.getEntry(i);
	}

	protected double f(int i)
	{
		return f.getEntry(i);
	}
	
	public static double sum(IndexedFunction f, int dimension, int exclude)
	{
		double sum = 0;
		for(int i = 0; i < dimension; ++i)
		{
			if(i != exclude)
			{
				sum += f.value(i);
			}
		}
		
		return sum;
	}
	
	public static double sum(IndexedFunction f, int dimension)
	{
		double sum = 0;
		for(int i = 0; i < dimension; ++i)
		{
			double v = f.value(i);
			sum += v;
		}
		
		return sum;
	}

	public static RealVector vectorByComponents(IndexedFunction f, int dimension)
	{
		RealVector grad = new ArrayRealVector(dimension);
		for(int i = 0; i < dimension; ++i)
			grad.setEntry(i, f.value(i));
		
		return grad; 
	}
	
	public static RealVector mask(RealVector v, IndexedBooleanFunction filter)
	{
		int dimension = v.getDimension();
		for (int i = 0; i < dimension; i++)
		{
			if(!filter.value(i))
				v.setEntry(i, 0);
		}
		return v;
	}
	
	
	//// Debug ////

	protected void debugHeader(AnalysisResults results)
	{
		Map<String, Integer> numbers = results.getBusNumbers();
		System.out.print("-,");
		
		int dimension = p.getDimension();
		debugLogNames(numbers, dimension, "p");
		debugLogNames(numbers, dimension, "q");
		debugLogNames(numbers, dimension, "e");
		debugLogNames(numbers, dimension, "f");
		debugLogNames(numbers, dimension, "l_p");
		debugLogNames(numbers, dimension, "l_q");
		debugLogNames(numbers, dimension, "gp");
		debugLogNames(numbers, dimension, "gq");
		debugLogNames(numbers, dimension, "gradL_p");
		debugLogNames(numbers, dimension, "gradL_q");
		debugLogNames(numbers, dimension, "gradL_e");
		debugLogNames(numbers, dimension, "gradL_f");
		System.out.print("C(pq),L(pql),t,c");
		System.out.println();
	}
	
	/**
	 * 
	 * @param k Iteration counter.
	 * @param t Step size.
	 */
	protected void debug(int k, double t)
	{
		System.out.print(k);
		System.out.print(',');
		int dimension = p.getDimension();
		
		// Power:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(p(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(q(i)));
			System.out.print(',');
		}
		
		// Voltage:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(e(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(f(i)));
			System.out.print(',');
		}
		
		// Lambda:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(lp(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(lq.getEntry(i)));
			System.out.print(',');
		}
		
		// Calculated Power Error:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(/*lp(i)**/gp(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(/*lq(i)**/gq(i)));
			System.out.print(',');
		}
		
		// Power derivatives:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(gradL_p(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(gradL_q(i)));
			System.out.print(',');
		}
		
		// Voltage derivatives:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(gradL_e(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(gradL_f(i)));
			System.out.print(',');
		}
		
		// Slack power:
//		System.out.print(format(pSlack()));
//		System.out.print(',');
//		System.out.print(format(qSlack()));
//		System.out.print(',');
		
		// Costs:
		System.out.print(format(cost()));
		System.out.print(',');
		System.out.print(format(L()));
		System.out.print(',');

		// Step size:
		System.out.print(t);//format(t));
		System.out.print(',');

		// Augmentation scale:
		System.out.print(c);
		System.out.print(',');
		
		System.out.println();
	}
	
	/**
	 * Lagrangian.
	 * @return
	 */
	protected double L()
	{
		int dimension = p.getDimension();
		
		return 
				cost()
				+sum(i -> lp(i)*gp(i), dimension)
				+sum(i -> lq(i)*gq(i), dimension)
				+0.5*c*sum(i -> gp(i)*gp(i), dimension, slackIndex)
				+0.5*c*sum(i -> gq(i)*gq(i), dimension, slackIndex);
	}

	private double cost()
	{
		int dimension = p.getDimension();
		return 
				sum(i -> generatorMask.getEntry(i) == 1 ? cost(i) : 0, dimension);
//				sum(i -> (1-e(i))*(1-e(i)) + f(i)*f(i), dimension);
	}

	protected double cost(int i)
	{
		return cost_i(p(i), q(i), p_max.getEntry(i));
	}

	protected double cost_i(double p_i, double q_i, double p_max_i)
	{
		double d = p_i - p_max_i;
		return 
				0.5*d*d
				+0.5*q_i*q_i;
	}
	
	private double cost(RealVector p, RealVector q)
	{
		int dimension = p.getDimension();
		return 
				sum(i -> generatorMask.getEntry(i) == 1 ? cost_i(p.getEntry(i), q.getEntry(i), p_max.getEntry(i)) : 0, dimension);
//				sum(i -> (1-e(i))*(1-e(i)) + f(i)*f(i), dimension);
	}

//	private static NumberFormat formatter = new DecimalFormat("0.######E0");
	public static String format(double d)
	{
		if(d < 1e6)
			return String.format("%.4f", d);
		else
			return String.format("%.4e", d);//formatter.format(d);
	}
	
	public static void logBusses(AnalysisResults results)
	{
		System.out.println("index,bus,e,f,p,q");
		Map<String, Integer> numbers = results.getBusNumbers();
		for (String name : numbers.keySet())
		{
			int index = numbers.get(name);
			System.out.print(index);
			System.out.print(',');
			System.out.print(name);
			System.out.print(',');
			Complex v = results.getBusVoltage(name);
			System.out.print(format(v.getReal()));
			System.out.print(',');
			System.out.print(format(v.getImaginary()));
			System.out.print(',');
			Complex s = results.getBusPower(name);
			System.out.print(format(s.getReal()));
			System.out.print(',');
			System.out.print(format(s.getImaginary()));
			System.out.println();
		}
	}
}