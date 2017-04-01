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
import static java.lang.Math.max;
import static java.lang.Math.min;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.Util;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Source;
import ellipsis.util.MatrixHelper;
import ellipsis.util.TeeOutputStream;

/**
 * Asynchronous version of {@link Sandbox014}.
 * @author bmillar
 *
 */
public class Sandbox015
{
	public static void main(String[] args)
	{
		new Sandbox015().run();
	}
	
	public static PrintStream out;
	static
	{
		try
		{
			out = new PrintStream(new TeeOutputStream(new FileOutputStream("/tmp/Sandbox015.csv"), System.out));
		} 
		catch (FileNotFoundException e)
		{
			throw new RuntimeException(e);
		}
	}
	
	
	//// Lambda expression interfaces ////
	
	public static interface IndexedFunction
	{
		double value(int i);
	}
	
	public static interface NumberFunction
	{
		double value(double d);
	}
	
	public static interface IndexedBooleanFunction
	{
		boolean value(int i);
	}
	
	public static interface NeighbourFunction
	{
		double value(Agent neighbour);
	}
	
	
	//// Global variables ////
	
	private static final double V_MAX = 1.05;
	private static final double V_MIN = 0.95;
	private static final double V_ARG_MIN = -Math.PI/4;
	private static final double V_ARG_MAX = Math.PI/4;
	private static final double CON_ERROR_TOLERANCE = 1e-3; // Maximum value of g(x) that's considered close enough to zero.
	private static final double MIN_CONSTRAINT_IMPROVEMENT = 0.9; // Minimum improvement percentage of g(x) constraint
	private static final double MAX_AUG_SCALE = 1e-1;
	private static final double AUG_SCALE_STEP = 1.1;
	private static final double X_STEPS_PER_ITERATION = 10;
	private static final int K = 1000;
	private static final double MAX_RAND_STEP_P = 1;
	private static final double MAX_RAND_STEP_Q = 0.5;
	private static final double MAX_RAND_STEP_E = 0.1;
	private static final double MAX_RAND_STEP_F = 0.05;
	private static final double MIN_STEP_SIZE = 0;//1e-9;
	
	// Switches (not final to avoid warnings):
	private static boolean ZERO_SMALL_G = false;
	private static boolean PROJECT_X = true;
	private static boolean USE_NEWTON = false;
	private static boolean FORMAT_NUMBERS = false;
	private static boolean RANDOMISE = false;
	private static boolean RANDOMISE_STEP_SIZE = false;
	
	protected RealVector p, q, e, f, lp, lq;
	protected RealMatrix G, B;
	protected RealVector p_max;
	private int slackIndex;
	protected RealVector generatorMask;
	protected Set<Agent> agents;
	protected Grid grid;
	protected double c = 1e-3;
	protected double cScale = 1;
	protected double t = 0;
	
	protected void run()
	{
		AnalysisResults results = init();
		loop();
		testResults(results);
		finish();
	}
	
	private void finish()
	{
		out.println();
		out.println(SimpleDateFormat.getDateTimeInstance().format((Calendar.getInstance().getTime())));
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
		lfa.setIterations(100);
		lfa.setTargetError(1e-6);
        AnalysisResults results2 = lfa.analyse();//complexVoltages()); // FIXME Passing in voltages is changing the load powers somehow in the results
        
        out.println("\nFinal Values:");
        out.println("index,bus,e,f,p,q,,|v|,|s|");
        for (String busName : busNumbers.keySet())
		{
			int i = busNumbers.get(busName);
        	double f_i = f(i);
			double e_i = e(i);
			String v_abs = format(Math.sqrt(e_i*e_i+f_i*f_i));
			double p_i = p(i);
			double q_i = q(i);
			String s_abs = format(Math.sqrt(p_i*p_i+q_i*q_i));
			out.println(i+","+busName+','+format(e_i)+','+format(f_i)+','+format(p_i)+','+format(q_i)+",,"+v_abs+','+s_abs);
        }
        out.println("\nAnalysis Results:");
//        GridDisplay.showInFrame(grid, grid.getBasePower(), grid.getBaseVoltage()).setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        logBusses(results2);
	}

//	private Complex[] complexVoltages()
//	{
//		int dimension = e.getDimension();
//		Complex[] v = new Complex[dimension];
//		for(int i = 0; i < dimension; ++i)
//		{
//			v[i] = new Complex(e(i), f(i));
//		}
//		return v;
//	}
	
	
	//// Agents ////
	
	public /*static*/ class Agent // FIXME
	{	
		public double p, q, e, f; // primal variables
		double p_new, q_new, e_new, f_new;
		public double p_max;
		public double lp, lq; // dual variables
		double lp_new, lq_new;
		public double gp, gq; // derived values
		protected double c;// = 1e-6;
		public boolean isSlack;
		public boolean isGenerator;
		public int index;
		private Map<Agent, Complex> neighbours = new HashMap<>();
		public double G, B; // Self admittances.
		public double lastT = 1;
		
		public Agent(int i)
		{
			this.index = i;
		}
		
		public void addNeighbour(Agent neighbour, Complex admittance)
		{
			neighbours.put(neighbour, admittance);
		}

		private boolean debug_backtrack = false;
//		private boolean debug_augscale = false;
		
		public void stepPrimal(int k)
		{
			if(isSlack)
				return;
	
if(k >= 9999 && index == 2)
{
	Util.nullop();
//	debug_backtrack = true;
}
			
			// Remember current constraint values:
			double gp_old = gp;//(p, q, e, f);
			double gq_old = gq;//(p, q, e, f);
			
			for(int i = 0; i < X_STEPS_PER_ITERATION; ++i)
				stepX();
			
			stepLambda();
			stepAugScale(k, gp_old, gq_old);
		}

		Random rand = new Random(0);
		int stepCount = 0;
		protected void stepX()
		{
			// Find x step direction:
			double dL_p = isGenerator ? gradL_p() : 0;
			double dL_q = isGenerator ? gradL_q() : 0;
			double dL_e = isSlack ? 0 : gradL_e();
			double dL_f = isSlack ? 0 : gradL_f();
			RealMatrix Hx = USE_NEWTON ? hessian() : MatrixUtils.createRealIdentityMatrix(4);
			RealVector x = new ArrayRealVector(new double[]{p, q, e, f});
			RealVector dL_x = new ArrayRealVector(new double[]{dL_p, dL_q, dL_e, dL_f});
			
			if(debug_backtrack)
			{
				compareGlobalLocal(dL_p, dL_q, dL_e, dL_f);
				System.out.println();
				checkGradient(dL_p, dL_q, dL_e, dL_f);
				System.out.println();
				plotImprovements();
				System.out.println();
				plotComparison();
			}
			
			// Backtrack to find best step size:
			double t = 1; // Step size - note that from Sandbox013 the smallest value is t=0.0313, and the average is 0.108 with c=0.01 I think.
			double alpha = 0.5; // Step update size.
			double beta = 0.5;
			double L_improvement;
			RealVector dx = MatrixHelper.solve(Hx, dL_x);// H_inv.operate(dL_x)
			if(dx == null)
				dx = dL_x; // Revert to gradient decent.
			dx = dx.mapMultiply(-1); // delta x = -H^-1 * grad_x L(x,lambda)
			double step = 
					alpha*dL_x.dotProduct(dx);
//					-alpha*(dL_p*dL_p + dL_q*dL_q + dL_e*dL_e + dL_f*dL_f);
			double minImprovement;

//			double dp = -dL_p;
//			double dq = -dL_q;
//			double de = -dL_e;
//			double df = -dL_f;
			double dp = dx.getEntry(0);
			double dq = dx.getEntry(1);
			double de = dx.getEntry(2);
			double df = dx.getEntry(3);
			
			do
			{
				L_improvement = lagrangianImprovement(t, dp, dq, de, df);
				minImprovement = t*step;

				if(debug_backtrack)
				{
					// Check improvement by L(x_t*dx)-L(x)
//					compareImprovement(dp, dq, de, df, t, L_improvement);
					double lowerBound = t*step/alpha;
					System.out.println(t+","+L_improvement+","+minImprovement+","+lowerBound);
				}
				
				t = t*beta;
			} while(L_improvement > minImprovement);// || debug_backtrack);
			
			// Randomise step size:
			if(RANDOMISE_STEP_SIZE)
			{
				t = Math.max(1.0, t + rand.nextDouble()*0.05);
			}
			
			// Limit step size:
			if(t < MIN_STEP_SIZE)
				t = MIN_STEP_SIZE;
			
			// Update state:
			RealVector x_new = x.add(dx.mapMultiply(t));
			p_new = x_new.getEntry(0);
			q_new = x_new.getEntry(1);
			e_new = x_new.getEntry(2);
			f_new = x_new.getEntry(3);
//			p_new = p - dL_p*t;
//			q_new = q - dL_q*t;
//			e_new = e - dL_e*t;
//			f_new = f - dL_f*t;
			
			if(debug_backtrack)
			{
				System.out.println("0,0,0,0");
				double _t = -2e-30;
				L_improvement = lagrangianImprovement(
					_t,
					dp, 
					dq, 
					de,
					df);
				double lowerBound = _t*step/alpha;
				System.out.println(_t+","+L_improvement+","+(_t*step)+","+lowerBound);
			}
			
			// Randomise (experimental):
			if(RANDOMISE)
			{
				p_new += MAX_RAND_STEP_P*(1-2*rand.nextDouble())*(10.0/(stepCount+10));
				q_new += MAX_RAND_STEP_Q*(1-2*rand.nextDouble())*(10.0/(stepCount+10));
				e_new += MAX_RAND_STEP_E*(1-2*rand.nextDouble())*(10.0/(stepCount+10));
				f_new += MAX_RAND_STEP_F*(1-2*rand.nextDouble())*(10.0/(stepCount+10));
			}
			
			// Project onto constrained set X:
			if(PROJECT_X)
			{
				Complex ef = projectVoltage(e_new, f_new);
				e_new = ef.getReal();
				f_new = ef.getImaginary();
				if(isGenerator)
				{
					Complex pq = projectPowers(p_new, q_new);
					p_new = pq.getReal();
					q_new = pq.getImaginary();
				}
			}
			
			// Update x:
			p = isGenerator ? p_new : p;
			q = isGenerator ? q_new : q;
			e = e_new;
			f = f_new;
			
			// Update cost function values:
			gp = gp();
			gq = gq();
			
			t /= beta; // Undo last increment.
			
			lastT = t;
			
			++stepCount;
		}

		private static final int P = 0;
		private static final int Q = 1;
		private static final int E = 2;
		private static final int F = 3;
		private RealMatrix hessian()
		{
			RealMatrix H = new Array2DRowRealMatrix(4, 4);
			H.setEntry(P, P, grad2L_pp());
			H.setEntry(P, Q, grad2L_pq());
			H.setEntry(P, E, grad2L_pe());
			H.setEntry(P, F, grad2L_pf());

			H.setEntry(Q, P, grad2L_qp());
			H.setEntry(Q, Q, grad2L_qq());
			H.setEntry(Q, E, grad2L_qe());
			H.setEntry(Q, F, grad2L_qf());

			H.setEntry(E, P, grad2L_ep());
			H.setEntry(E, Q, grad2L_eq());
			H.setEntry(E, E, grad2L_ee());
			H.setEntry(E, F, grad2L_ef());

			H.setEntry(F, P, grad2L_fp());
			H.setEntry(F, Q, grad2L_fq());
			H.setEntry(F, E, grad2L_fe());
			H.setEntry(F, F, grad2L_ff());
			
			return H;
		}

		private double grad2L_pp()
		{
			return 1 + c;
		}

		private double grad2L_pq()
		{
			return 0;
		}

		private double grad2L_pe()
		{
			return -c*gradGp_e();
		}

		private double grad2L_pf()
		{
			return -c*gradGp_f();
		}

		private double grad2L_qp()
		{
			return 0;
		}

		private double grad2L_qq()
		{
			return 1 + c;
		}

		private double grad2L_qe()
		{
			return -c*gradGq_e();
		}

		private double grad2L_qf()
		{
			return -c*gradGq_f();
		}

		private double grad2L_ep()
		{
			return -c*sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) + 2*G*e;
		}

		private double grad2L_eq()
		{
			return -c*sum(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent)) - 2*B*e;
		}

		private double grad2L_ee()
		{
			return
					2*lp*G - 2*lq*B
					+c*sum(neighbours, agent -> 
							  gradGp_e(agent)*(agent.e*G(agent) + agent.f*B(agent))
							+ gradGq_e(agent)*(agent.f*G(agent) - agent.e*B(agent))
						)
					+c*gradGp_e()*sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent))
					+c*(gradGp_e()*2*G*e + 2*G*gp())
					+c*gradGq_e()*sum(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent))
					-c*(gradGq_e()*2*B*e + 2*B*gq())
					;
		}

		private double gradGp_e(Agent agent)
		{
			return agent.e*G(agent) + agent.f*B(agent);
		}

		private double gradGp_f(Agent agent)
		{
			return agent.f*G(agent) - agent.e*B(agent);
		}

		private double gradGq_e(Agent agent)
		{
			return agent.f*G(agent) - agent.e*B(agent);
		}

		private double gradGq_f(Agent agent)
		{
			return -agent.e*G(agent) - agent.f*B(agent);
		}

		private double grad2L_ef()
		{
			return
					c*sum(neighbours, agent -> gradGp_f(agent)*(agent.e*G(agent) + agent.f*B(agent)) + gradGq_f(agent)*(agent.f*G(agent) - agent.e*B(agent)))
					+c*gradGp_f()*(sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) + 2*G*e)
					+c*gradGq_f()*(sum(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent)) - 2*B*e);
		}

		private double grad2L_fp()
		{
			return -c*sum(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent)) + 2*G*f;
		}

		private double grad2L_fq()
		{
			return -c*sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) - 2*B*f;
		}

		private double grad2L_fe()
		{
			return
					c*sum(neighbours, agent -> gradGp_e(agent)*(agent.f*G(agent) - agent.e*B(agent)) + gradGq_e(agent)*(-agent.e*G(agent) - agent.f*B(agent)))
					+c*gradGp_e()*(sum(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent)) + 2*G*f)
					+c*gradGq_e()*(sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) - 2*B*f);
		}

		private double grad2L_ff()
		{
			return
					2*lp*G - 2*lq*B
					+c*sum(neighbours, agent -> 
							  gradGp_f(agent)*(agent.f*G(agent) - agent.e*B(agent))
							+ gradGq_f(agent)*(-agent.e*G(agent) - agent.f*B(agent))
						)
					+c*gradGp_f()*sum(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent))
					+c*(gradGp_f()*2*G*f + 2*G*gp())
					+c*gradGq_f()*sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent))
					-c*(gradGq_f()*2*B*f + 2*B*gq())
					;
		}
		
		private double gradGp_e()
		{
			return sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) + 2*G*e;
		}
		
		private double gradGp_f()
		{
			return sum(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent)) + 2*G*f;
		}
		
		private double gradGq_e()
		{
			return sum(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent)) - 2*B*e;
		}
		
		private double gradGq_f()
		{
			return sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) - 2*B*f;
		}

		protected void stepLambda()
		{
			// Maximise (l := l - c*g(x), ref http://en.wikipedia.org/wiki/Augmented_Lagrangian_method
			//           and Bertsekas' book):
			double dlp = isSlack ? 0 : gp;
			lp = lp + dlp*c;
			double dlq = isSlack ? 0 : gq;
			lq = lq + dlq*c;
			
			lp = project(lp, -1e14, 1e14);
			lq = project(lq, -1e14, 1e14);
		}

		protected void stepAugScale(int k, double gp_old, double gq_old)
		{
			// Update augmentation scale:
			if(Math.abs(gp) > MIN_CONSTRAINT_IMPROVEMENT*Math.abs(gp_old) || Math.abs(gq) > MIN_CONSTRAINT_IMPROVEMENT*Math.abs(gq_old))
			{
				if(c < MAX_AUG_SCALE)
					c *= AUG_SCALE_STEP;
			}
			
			// experimental
//			if(Math.abs(gp) > MIN_CONSTRAINT_IMPROVEMENT*Math.abs(gp_old) || Math.abs(gq) > MIN_CONSTRAINT_IMPROVEMENT*Math.abs(gq_old))
//			{
//				if(cScale < MAX_AUG_SCALE)
//					cScale *= AUG_SCALE_STEP;
//			}
//			c = cScale*(1+Math.sin(2*Math.PI*(k/10.0)));
			
//			if(debug_augscale)
//			{
//				out.println(gp(p, q, e, f).getNorm()+" ?> "+delta+"*"+gp.getNorm()+"="+delta*gp.getNorm());
//				out.println(gq(p, q, e, f).getNorm()+" ?> "+delta+"*"+gq.getNorm()+"="+delta*gq.getNorm());
//			}
		}

		private Complex projectVoltage(double e, double f)
		{
			if(e < 0)
				e = 0;
			
			// Clamp |v|:
			Complex v = new Complex(e, f);
			double abs = v.abs();
			if(abs > V_MAX)
			{
				v = v.multiply(V_MAX/abs);
			}
			else if(abs < V_MIN)
			{
				v = v.multiply(V_MIN/abs);
			}
			
			// Clamp angle:
			abs = v.abs();
			double arg = v.getArgument();
			if(arg < V_ARG_MIN)
			{
				v = ComplexUtils.polar2Complex(abs, V_ARG_MIN);
			}
			else if(arg > V_ARG_MAX)
			{
				v = ComplexUtils.polar2Complex(abs, V_ARG_MAX);
			}
			
			return v;
		}

		private Complex projectPowers(double p, double q)
		{
			if(p < 0)
				p = 0;
			
			Complex s = new Complex(p, q);
			double abs = s.abs();
			double p_max_i = p_max;
			if(abs > p_max_i)
			{
				s = s.multiply(p_max_i/abs);
			}
			
			return s;
		}

		/**
		 * This is equivalent to L(X+tdx) - L(x) but is calculable within the agent's neighbourhood.
		 * @param t
		 * @param dp
		 * @param dq
		 * @param de
		 * @param df
		 * @return
		 */
		private double lagrangianImprovement(double t, double dp, double dq, double de, double df)
		{
			double dc = cost(p+t*dp, q+t*dq) - cost(p, q);
//			double gp_dx = gp(t*dp, t*dq, t*de, t*df);
//			double gq_dx = gq(t*dp, t*dq, t*de, t*df);
//			double dgp = gp_dx + 2*t*G*(e*de + f*df);
//			double dgq = gq_dx - 2*t*B*(e*de + f*df);
			
			
			double constraintTerms = lagrangianImprovement_constraintTerms(t, dp, dq, de, df);
			double augTermP = lagrangianImprovement_augTermP(t, dp, dq, de, df);
			double augTermQ = lagrangianImprovement_augTermQ(t, dp, dq, de, df);
			double augTerm = augTermP + augTermQ;
			return 
//					dc
//					+lp*dgp
//					+lq*dgq
//					+0.5*c*(dgp*dgp + 2*gp*dgp)
//					+0.5*c*(dgq*dgq + 2*gq*dgq);
					dc
					+constraintTerms
					+augTerm;
		}

		public double lagrangianImprovement_augTermP(double t, double dp, double dq, double de, double df)
		{
			return 0.5*c*sum(neighbours, agent -> agent.isSlack ? 0 : dgp(agent, t, dp, dq, de, df)*(dgp(agent, t, dp, dq, de, df) + 2*agent.gp()))
								+0.5*c*( dgpSelf(t, dp, dq, de, df)*(dgpSelf(t, dp, dq, de, df) + 2*gp()) );
		}

		public double lagrangianImprovement_augTermQ(double t, double dp, double dq, double de, double df)
		{
			return 
					 0.5*c*sum(neighbours, agent -> agent.isSlack ? 0 : dgq(agent, t, dp, dq, de, df)*(dgq(agent, t, dp, dq, de, df) + 2*agent.gq()))
					+0.5*c*(dgqSelf(t, dp, dq, de, df)*(dgqSelf(t, dp, dq, de, df) + 2*gq()));
		}

		public double lagrangianImprovement_constraintTerms(double t,
				double dp, double dq, double de, double df)
		{
			return 
					sum(neighbours, agent -> agent.isSlack ? 0 : agent.lp*dgp(agent, t, dp, dq, de, df))
					+lp*dgpSelf(t, dp, dq, de, df)
					+sum(neighbours, agent -> agent.isSlack ? 0 : agent.lq*dgq(agent, t, dp, dq, de, df))
					+lq*dgqSelf(t, dp, dq, de, df);
		}

		private double dgpSelf(double t, double dp, double dq, double de, double df)
		{
			return
					t*t*(de*de + df*df)*G - t*dp
					+t*(e*de*G + f*df*G + f*de*B - e*df*B)
					+t*sum(neighbours, agent -> de*agent.e*G(agent) + df*agent.f*G(agent) + df*agent.e*B(agent) - de*agent.f*B(agent))
					+t*(de*e*G + df*f*G + df*e*B - de*f*B); // This last term is needed since neighbours does not include this agent
		}

		private double dgqSelf(double t, double dp, double dq, double de, double df)
		{
			return
					-t*t*(de*de + df*df)*B - t*dq
					+t*(f*de*G - e*df*G - e*de*B - f*df*B)
					+t*sum(neighbours, agent -> df*agent.e*G(agent) - de*agent.f*G(agent) - de*agent.e*B(agent) - df*agent.f*B(agent))
					+t*(df*e*G - de*f*G - de*e*B - df*f*B); // This last term is needed since neighbours does not include this agent
		}

		private double dgp(Agent agent, double t, double dp, double dq, double de, double df)
		{
			double G_in = G(agent);
			double B_in = B(agent);
			return t*(agent.e*de*G_in + agent.f*df*G_in + agent.f*de*B_in - agent.e*df*B_in);
		}

		private double dgq(Agent agent, double t, double dp, double dq, double de, double df)
		{
			double G_in = G(agent);
			double B_in = B(agent);
			return t*(agent.f*de*G_in - agent.e*df*G_in - agent.e*de*B_in - agent.f*df*B_in);
		}

		private double project(double v, double low, double high)
		{
			return min(max(low, v), high);
		}

		private double gp()
		{
			double d = gp(p, q, e, f);
			
			return ZERO_SMALL_G && (Math.abs(d) < CON_ERROR_TOLERANCE) ? 0 : d;
		}
		
		/**
		 * WARNING: This still uses neighbour values, only this agent's values are passed in.
		 * @param p
		 * @param q
		 * @param e
		 * @param f
		 * @return
		 */
		private double gp(double p, double q, double e, double f)
		{
			double sum = 
					0;
					//sum(neighbours, agent -> e*agent.e*G(agent) + f*agent.f*G(agent) + f*agent.e*B(agent) - e*agent.f*B(agent));
			for (Agent agent : neighbours.keySet())
			{
				sum += e*agent.e*G(agent) + f*agent.f*G(agent) + f*agent.e*B(agent) - e*agent.f*B(agent);
			}
			sum += e*e*G + f*f*G; // self
			return 	sum	- p;
		}
		
		private double gq()
		{
			double d = gq(p, q, e, f);
			return ZERO_SMALL_G && (Math.abs(d) < CON_ERROR_TOLERANCE) ? 0 : d;
		}
		
		private double gq(double p, double q, double e, double f)
		{
			double sum1 = sum(neighbours, agent -> f*agent.e*G(agent) - e*agent.f*G(agent) - e*agent.e*B(agent) - f*agent.f*B(agent));
			double sum = sum1
					- e*e*B - f*f*B;
			return 
				sum
				- q;
		}
		
		/**
		 * c(x) = 0.5*(p - p_max)^2 + 0.5*q^2
		 * @param p Real power output.
		 * @param q Reactive power output.
		 * @return
		 */
		private double cost(double p, double q)
		{
			double d = p-p_max;
			return 0.5*d*d + 0.5*q*q;
		}

		protected double gradL_p()
		{
			return 	gradC_p() - lp

					// Penalty:
					-c*gp()
					;
		}
		
		/**
		 * Assuming C_i(p_i) = 0.5*(p_i - p_max)^2
		 * => dC_i/dp_i = p_i - p_max
		 * @return
		 */
		private double gradC_p()
		{
			return p - p_max;
		}

		protected double gradL_q()
		{
			return gradC_q() - lq

					// Penalty:
					-c*gq() // FIXME
					;
		}
		
		/**
		 * Assuming C_i(q_i) = 0.5*q_i^2
		 * => dC_i/dp_i = q_i
		 * @param i
		 * @return
		 */
		private double gradC_q()
		{
			return q;
		}

		protected double gradL_e()
		{
			double constraintTerms = 
						 sum(neighbours, agent -> agent.lp*(agent.e*G(agent) + agent.f*B(agent)))
						+sum(neighbours, agent -> agent.lq*(agent.f*G(agent) - agent.e*B(agent)))
						+lp*(sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent))
						     +2*G*e)
						+lq*(sum(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent))
						     -2*B*e);

			double pt1 = c*sum(neighbours, agent -> agent.isSlack ? 0 : 
								 agent.gp()*(agent.e*G(agent) + agent.f*B(agent))
								+agent.gq()*(agent.f*G(agent) - agent.e*B(agent)));
			double pt2 = c*gp()*(sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) 
					           +2*e*G);
			double pt3 = c*gq()*(sum(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent)) 
					           -2*e*B);
			double penaltyTerms = pt1
								+
								(isSlack ? 0 : 
								(
									+pt2
									+pt3
								));
			
			return
				constraintTerms
				+penaltyTerms
				;
		}

		protected double gradL_f()
		{
			double constraintTerms = 
						 sum(neighbours, agent -> agent.lp*(agent.f*G(agent) - agent.e*B(agent)))
						+sum(neighbours, agent -> agent.lq*(-agent.e*G(agent) - agent.f*B(agent)))
						+lp*(sum(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent))
						     +2*G*f)
						+lq*(sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent))
						     -2*B*f);

			double pt1 = c*sum(neighbours, agent -> agent.isSlack ? 0 : 
							   agent.gp()*(agent.f*G(agent) - agent.e*B(agent))
							   +agent.gq()*(-agent.e*G(agent) - agent.f*B(agent)));
			double pt2 = c*gp()*(sum(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent)) 
				               +2*f*G);
			double pt3 = c*gq()*(sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) 
				               -2*f*B);
			double penaltyTerms = pt1
								+
								(isSlack ? 0 : 
								(
									+pt2
									+pt3
								));
			
//out.print(c*gq()+",*(,");
//for(Agent agent : neighbours.keySet())
//{
//	out.print(agent.e*G(agent) - agent.f*B(agent));
////	out.println(
////	 agent.gp()+",*( +,"+(agent.f*G(agent) - agent.e*B(agent))+",),"+
////	+agent.gq()+",*( +,"+(-agent.e*G(agent) - agent.f*B(agent))+",)");
//}
//out.println(",) +,"+(-2*f*B));
			
			return
				constraintTerms
				+penaltyTerms
				;
		}

		private double G(Agent agent)
		{
			double G_ij = neighbours.get(agent).getReal();
			return G_ij;
		}

		private double B(Agent agent)
		{
			return neighbours.get(agent).getImaginary();
		}

		public void init()
		{
			gp = gp();
			gq = gq();
			
			// Find unconstrained optimum of x as the starting point:
			if(isGenerator)
			{
				p = p_max;
				q = 0;
			}
		}
		
		public String toString()
		{
			StringBuffer buf = new StringBuffer();
			buf.append("Agent ");
			buf.append(index);
			buf.append(": \n");
			buf.append("\tp=");
			buf.append(p);
			buf.append("\n\tq=");
			buf.append(q);
			buf.append("\n\te=");
			buf.append(e);
			buf.append("\n\tf=");
			buf.append(f);
			buf.append("\n\tlp=");
			buf.append(lp);
			buf.append("\n\tlq=");
			buf.append(lq);
			buf.append("\n\tc=");
			buf.append(c);
			buf.append("\n\tt=");
			buf.append(lastT);
			buf.append("\n");
				
			return buf.toString();
		}
		
		
		//// Debug ////
		
		public Object showNeighbours = new Object()
		{
			public String toString()
			{
				StringBuffer buf = new StringBuffer();

				for (Agent n : neighbours.keySet())
				{
					buf.append(n);
					buf.append("\tG_ij=");
					buf.append(neighbours.get(n).getReal());
					buf.append("\n\tB_ij=");
					buf.append(neighbours.get(n).getImaginary());
					buf.append("\n\n");
				}
				
				return buf.toString();
			}
		};

		protected void compareGlobalLocal(double dp, double dq, double de,
				double df)
		{
			System.out.println("Compare Global/Local:");
			
			double globalDp = Sandbox015.this.gradL_p(index);
			System.out.println(globalDp + ",?=," + dp);

			double globalDq = Sandbox015.this.gradL_q(index);
			System.out.println(globalDq + ",?=," + dq);
			
			double globalDe = Sandbox015.this.gradL_e(index);
			System.out.println(globalDe + ",?=," + de);

			double globalDf = Sandbox015.this.gradL_f(index);
			System.out.println(globalDf + ",?=," + df);
		}

		protected void compareImprovement(double dp, double dq, double de,
				double df, double t, double L_improvement)
		{
			RealVector _p = new ArrayRealVector(Sandbox015.this.p);
			RealVector _q = new ArrayRealVector(Sandbox015.this.q);
			RealVector _e = new ArrayRealVector(Sandbox015.this.e);
			RealVector _f = new ArrayRealVector(Sandbox015.this.f);
			double Lx = lagrangian(_p, _q, _e, _f, Sandbox015.this.lp, Sandbox015.this.lq);
			_p.setEntry(index, p-t*dp);
			_q.setEntry(index, q-t*dq);
			_e.setEntry(index, e-t*de);
			_f.setEntry(index, f-t*df);
			double Ldx = lagrangian(_p, _q, _e, _f, Sandbox015.this.lp, Sandbox015.this.lq);
			System.out.println("L(x_t*dx)-L(x) =, "+(Ldx-Lx)+", L_imp =, "+L_improvement);
		}

		private void plotComparison()
		{
			System.out.println(
				    "Comparison:\ndx,imp_g(x)_e,imp_g(x)_f,g(x+de)-g(x),g(x_df)-g(x)"
					+ ",imp_augP(x)_e,imp_augP(x)_f,augP(x+de)-augP(x),augP(x_df)-augP(x)"
					+ ",imp_augQ(x)_e,imp_augQ(x)_f,augQ(x+de)-augQ(x),augQ(x_df)-augQ(x)"
					+ ",imp_sum(x)_e,imp_sum(x)_f,L_sum(x+de)-L_sum(x),L_sum(x_df)-L_sum(x)"
					+ ",imp(x)_e,imp(x)_f,L(x+de)-L(x),L(x_df)-L(x)");
			for(double d = -1e-1; d <= 1.0000001e-1; d += 1e-2)
			{
				// Local calcs:
				double _de_cons = lagrangianImprovement_constraintTerms(1, 0, 0, -d, 0);
				double _df_cons = lagrangianImprovement_constraintTerms(1, 0, 0, 0, -d);
				
				double _de_augP = lagrangianImprovement_augTermP(1, 0, 0, -d, 0);
				double _df_augP = lagrangianImprovement_augTermP(1, 0, 0, 0, -d);
				
				double _de_augQ = lagrangianImprovement_augTermQ(1, 0, 0, -d, 0);
				double _df_augQ = lagrangianImprovement_augTermQ(1, 0, 0, 0, -d);

				double _de = lagrangianImprovement(1, 0, 0, -d, 0);
				double _df = lagrangianImprovement(1, 0, 0, 0, -d);
				
				// Global calcs:
				RealVector _p = new ArrayRealVector(Sandbox015.this.p);
				RealVector _q = new ArrayRealVector(Sandbox015.this.q);
				RealVector _e = new ArrayRealVector(Sandbox015.this.e);
				RealVector _f = new ArrayRealVector(Sandbox015.this.f);
				RealVector gp_global = Sandbox015.this.gp(_p, _q, _e, _f);
				RealVector gq_global = Sandbox015.this.gq(_p, _q, _e, _f);
				gp_global.setEntry(slackIndex, 0);
				gq_global.setEntry(slackIndex, 0);
				double Lx_augP = 0.5*c*gp_global.dotProduct(gp_global);
				double Lx_augQ = 0.5*c*gq_global.dotProduct(gq_global);
				double Lx_lg = Sandbox015.this.lp.dotProduct(gp_global) + Sandbox015.this.lq.dotProduct(gq_global);
				double Lx = lagrangian(_p, _q, _e, _f, Sandbox015.this.lp, Sandbox015.this.lq);

				_e.setEntry(index, e-d);
				gp_global = Sandbox015.this.gp(_p, _q, _e, _f);
				gq_global = Sandbox015.this.gq(_p, _q, _e, _f);
				gp_global.setEntry(slackIndex, 0);
				gq_global.setEntry(slackIndex, 0);
				double Lx_augP_e = 0.5*c*gp_global.dotProduct(gp_global);
				double Lx_augQ_e = 0.5*c*gq_global.dotProduct(gq_global);
				double Lx_lg_e = Sandbox015.this.lp.dotProduct(gp_global) + Sandbox015.this.lq.dotProduct(gq_global);
				double Lx_e = lagrangian(_p, _q, _e, _f, Sandbox015.this.lp, Sandbox015.this.lq);
				
				_e.setEntry(index, e);
				_f.setEntry(index, f-d);
				gp_global = Sandbox015.this.gp(_p, _q, _e, _f);
				gq_global = Sandbox015.this.gq(_p, _q, _e, _f);
				gp_global.setEntry(slackIndex, 0);
				gq_global.setEntry(slackIndex, 0);
				double Lx_augP_f = 0.5*c*gp_global.dotProduct(gp_global);
				double Lx_augQ_f = 0.5*c*gq_global.dotProduct(gq_global);
				double Lx_lg_f = Sandbox015.this.lp.dotProduct(gp_global) + Sandbox015.this.lq.dotProduct(gq_global);
				double Lx_f = lagrangian(_p, _q, _e, _f, Sandbox015.this.lp, Sandbox015.this.lq);
				
				System.out.println(d+","+_de_cons+","+_df_cons+","+(Lx_lg_e-Lx_lg)+","+(Lx_lg_f-Lx_lg)
						+","+_de_augP+","+_df_augP+","+(Lx_augP_e-Lx_augP)+","+(Lx_augP_f-Lx_augP)
						+","+_de_augQ+","+_df_augQ+","+(Lx_augQ_e-Lx_augQ)+","+(Lx_augQ_f-Lx_augQ)
						+","+(_de_cons+_de_augP+_de_augQ)+","+(_df_cons+_df_augP+_df_augQ)+","
							+((Lx_lg_e+Lx_augP_e+Lx_augQ_e)-(Lx_lg+Lx_augP+Lx_augQ))+","+((Lx_lg_f+Lx_augP_f+Lx_augQ_f)-(Lx_lg+Lx_augP+Lx_augQ))
						+","+_de+","+_df+","+(Lx_e-Lx)+","+(Lx_f-Lx));
			}
		}

		private void plotImprovements()
		{
			System.out.println("Improvements:\ndx,imp_e,imp_f,L(x+de)-L(x),L(x_df)-L(x)");
//			for(double d = -1e-12; d < 1e-12; d += 1e-13)
			for(double d = -1e-14; d <= 1.0000001e-14; d += 1e-15)
			{
				double _de = lagrangianImprovement(1, 0, 0, -d, 0);
				double _df = lagrangianImprovement(1, 0, 0, 0, -d);
				
				RealVector _p = new ArrayRealVector(Sandbox015.this.p);
				RealVector _q = new ArrayRealVector(Sandbox015.this.q);
				RealVector _e = new ArrayRealVector(Sandbox015.this.e);
				RealVector _f = new ArrayRealVector(Sandbox015.this.f);
				double Lx = lagrangian(_p, _q, _e, _f, Sandbox015.this.lp, Sandbox015.this.lq);
				_e.setEntry(index, e-d);
				double Ldx_e = lagrangian(_p, _q, _e, _f, Sandbox015.this.lp, Sandbox015.this.lq);
				
				_e.setEntry(index, e);
				_f.setEntry(index, f-d);
				double Ldx_f = lagrangian(_p, _q, _e, _f, Sandbox015.this.lp, Sandbox015.this.lq);
				
				System.out.println(d+","+_de+","+_df+","+(Ldx_e-Lx)+","+(Ldx_f-Lx));
			}
			
//			for(double d = -1e-12; d < 1e-12; d += 1e-13)
//			{
//				double _imp = lagrangianImprovement(1, 0, 0, -d, -d);
//				RealVector _p = new ArrayRealVector(Sandbox015.this.p);
//				RealVector _q = new ArrayRealVector(Sandbox015.this.q);
//				RealVector _e = new ArrayRealVector(Sandbox015.this.e);
//				RealVector _f = new ArrayRealVector(Sandbox015.this.f);
//				double Lx = lagrangian(_p, _q, _e, _f, Sandbox015.this.lp, Sandbox015.this.lq);
//				_e.setEntry(index, e-d);
//				_f.setEntry(index, f-d);
//				double Ldx = lagrangian(_p, _q, _e, _f, Sandbox015.this.lp, Sandbox015.this.lq);
//				out.println(d+","+_imp+","+(Ldx-Lx));
//			}
		}

		private void checkGradient(double dp, double dq, double de, double df)
		{
			double delta = 1e-24;
			
			System.out.println("Check Gradient:\nindex,dp,,dq,,de,,df");
			
			// Estimate p slope:
			double left = lagrangianImprovement(1, -delta, 0, 0, 0);
			double right = lagrangianImprovement(1, delta, 0, 0, 0);
			double pSlope = (right-left)/(2*delta);

			// Estimate q slope:
			left = lagrangianImprovement(1, 0, -delta, 0, 0);
			right = lagrangianImprovement(1, 0, delta, 0, 0);
			double qSlope = (right-left)/(2*delta);
			
			// Estimate e slope:
			left = lagrangianImprovement(1, 0, 0, -delta, 0);
			right = lagrangianImprovement(1, 0, 0, delta, 0);
			double eSlope = (right-left)/(2*delta);

			// Estimate f slope:
			left = lagrangianImprovement(1, 0, 0, 0, -delta);
			right = lagrangianImprovement(1, 0, 0, 0, delta);
			double fSlope = (right-left)/(2*delta);
			
			// Log:
			System.out.println(index
					+","+dp+","+pSlope
					+","+dq+","+qSlope
					+","+de+","+eSlope
					+","+df+","+fSlope);
		}
	} // End Agent
	
	
	//// Setup ////

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
        
        // Create agents:
        initAgents();
        
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

	public static RealVector p0(AnalysisResults results)
	{
		Map<String, Integer> numbers = results.getBusNumbers();
		double[] p0 = new double[numbers.size()];
		
		for (String name : numbers.keySet())
		{
			int i = numbers.get(name);
			p0[i] = results.getBusPower(name).getReal();
		}
		
		return new ArrayRealVector(p0);
	}

	public static RealVector q0(AnalysisResults results)
	{
		Map<String, Integer> numbers = results.getBusNumbers();
		double[] q0 = new double[numbers.size()];
		
		for (String name : numbers.keySet())
		{
			int i = numbers.get(name);
			q0[i] = results.getBusPower(name).getImaginary();
		}
		
		return new ArrayRealVector(q0);
	}

	public static RealVector maxPower(AnalysisResults results, Grid grid)
	{
		Map<String, Integer> numbers = results.getBusNumbers();
		RealVector powers = new ArrayRealVector(numbers.size());
		for (String name : numbers.keySet())
		{
			int index = numbers.get(name);
			Bus bus = grid.getBus(name);
			Collection<Source> sources = bus.getSources();
			if(sources.size() > 0)
			{
				Source s = sources.iterator().next();
				powers.setEntry(index, s.getPmax()/grid.getBasePower());
			}
		}
		return powers;
	}
	
	protected void initAgents()
	{
		// Make agents:
		agents = new HashSet<>();
		int dimension = p.getDimension();
		for(int i = 0; i < dimension; ++i)
		{
			Agent agent = new Agent(i);
			agent.c = this.c;
			agent.p = p(i);
			agent.q = q(i);
			agent.e = e(i);
			agent.f = f(i);
			agent.p_max = p_max.getEntry(i);
			agent.lp = lp(i);
			agent.lq = lq(i);
			agent.isGenerator = generatorMask.getEntry(i) == 1;
			agent.isSlack = i == slackIndex;
			agents.add(agent);
		}
		
		// Setup Neighbours:
		for (Agent agent_i : agents)
		{
			int i = agent_i.index;
			for (Agent agent_j : agents)
			{
				int j = agent_j.index;
				Complex Y_ij = new Complex(G(i,j), B(i,j));
				if(!Y_ij.equals(Complex.ZERO))
				{
					if(i == j)
					{
						agent_i.G = Y_ij.getReal();
						agent_i.B = Y_ij.getImaginary();
					}
					else
					{
						agent_i.addNeighbour(agent_j, Y_ij);
					}
				}
			}
		}
		
		for (Agent agent : agents)
			agent.init();
	}
	
	
	//// Iterations ////

	public void loop()
	{
		for(int k = 0; k < K; ++k)
		{
			for (Agent agent : agents)
			{
				agent.stepPrimal(k);
				update();
			}
			
			update();

			// Log:
			debug(k, t);
		}
	}
	
	/**
	 * Update global variables from agent values.
	 */
	private void update()
	{
		double count = 0;
		t = 0;
		c = 0;
		for (Agent agent : agents)
		{
			int i = agent.index;
			p.setEntry(i, agent.p);
			q.setEntry(i, agent.q);
			e.setEntry(i, agent.e);
			f.setEntry(i, agent.f);
			lp.setEntry(i, agent.lp);
			lq.setEntry(i, agent.lq);

			if(!agent.isSlack)
			{
				c += agent.c;
				t += agent.lastT;
				++count;
			}
		}
		t /= count;
		c /= count;
	}
	
	
	//// Convenience methods ////

	public static String format(double d)
	{
		if(!FORMAT_NUMBERS)
			return ""+d;
		
		double abs = Math.abs(d);
		if(abs < 1e-24)
			return "0.0";
		else if(abs > 1e-3 && abs < 1e6)
			return String.format("%.4f", d);
		else
			return String.format("%.4e", d);
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
			sum += f.value(i);
		}
		
		return sum;
	}
	
	public static double sum(double[] ds, NumberFunction f)
	{
		double sum = 0;
		for(int i = 0; i < ds.length; ++i)
		{
			sum += f.value(ds[i]);
		}
		
		return sum;
	}
	
	public static double sum(Map<Agent, Complex> neighbours, NeighbourFunction f)
	{
		double sum = 0;
		for (Agent agent : neighbours.keySet())
		{
			double v = f.value(agent);
			sum += v;
		}
		
		return sum;
	}

	public static double norm(double... ds)
	{
		return Math.sqrt(sum(ds, d -> d*d));
	}

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
	
	
	//// Debug ////

	protected void debugHeader(AnalysisResults results)
	{
		Map<String, Integer> numbers = results.getBusNumbers();
		out.print("-,");
		
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
		out.print("C(pq),L(pql),t,c,");
		debugLogNames(numbers, dimension, "inf");
		out.println();
	}

	public static void debugLogNames(Map<String, Integer> numbers, int dimension, String type)
	{
		for(int i = 0; i < dimension; ++i)
		{
			for (String name : numbers.keySet())
			{
				int index = numbers.get(name);
				if(index == i)
				{
					out.print(name);
					out.print('(');
					out.print(type);
					out.print("),");
					break;
				}
			}
		}
	}
	
	/**
	 * 
	 * @param k Iteration counter.
	 * @param t 
	 */
	protected void debug(int k, double t)
	{
		out.print(k);
		out.print(',');
		int dimension = p.getDimension();
		
		// Power:
		for(int i = 0; i < dimension; ++i)
		{
			out.print(format(p(i)));
			out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			out.print(format(q(i)));
			out.print(',');
		}
		
		// Voltage:
		for(int i = 0; i < dimension; ++i)
		{
			out.print(format(e(i)));
			out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			out.print(format(f(i)));
			out.print(',');
		}
		
		// Lambda:
		for(int i = 0; i < dimension; ++i)
		{
			out.print(format(lp(i)));
			out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			out.print(format(lq.getEntry(i)));
			out.print(',');
		}
		
		// Calculated Power Error:
		for(int i = 0; i < dimension; ++i)
		{
			out.print(format(gp(i)));
			out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			out.print(format(gq(i)));
			out.print(',');
		}
		
		// Power derivatives:
		for(int i = 0; i < dimension; ++i)
		{
			out.print(format(gradL_p(i)));
			out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			out.print(format(gradL_q(i)));
			out.print(',');
		}
		
		// Voltage derivatives:
		for(int i = 0; i < dimension; ++i)
		{
			out.print(format(gradL_e(i)));
			out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			out.print(format(gradL_f(i)));
			out.print(',');
		}
		
		// Slack power:
//		out.print(format(pSlack()));
//		out.print(',');
//		out.print(format(qSlack()));
//		out.print(',');
		
		// Costs:
		out.print(format(cost()));
		out.print(',');
		out.print(format(L()));
		out.print(',');

		// Step size:
		out.print(t);
		out.print(',');

		// Augmentation scale:
		out.print(c);
		out.print(',');

		// Controllability:
		for(int i = 0; i < dimension; ++i)
		{
			Agent agent = findAgent(i);
			out.print(format(controllability(agent)));
			out.print(',');
		}
		out.print(',');
		
		out.println();
	}

	private Agent findAgent(int i)
	{
		for (Agent agent : agents)
		{
			if(agent.index == i)
				return agent;
		}
		return null;
	}

	private double controllability(Agent agent)
	{
		double gradCNorm = norm(agent.gradC_p(), agent.gradC_q());
		double gradGNorm = norm(agent.gradGp_e(), agent.gradGp_f(), agent.gradGq_e(), agent.gradGq_f());
		double gNorm = norm(agent.gp(), agent.gq());
		double lNorm = norm(agent.lp, agent.lq);
		return gradCNorm/(gradGNorm*(lNorm+agent.c*gNorm));
	}

	protected double gradL_p(int i)
	{
		return 	gradF_p(i) - lp(i)

				// Penalty:
				-c*gp(i)
				;
	}
	
	/**
	 * Assuming F_i(p_i) = 0.5*(p_i - p_max)^2
	 * => dF_i/dp_i = p_i - p_max
	 * @param i
	 * @return
	 */
	private double gradF_p(int i)
	{
		return p(i) - p_max.getEntry(i);
	}

	protected double gradL_q(int i)
	{
		return gradF_q(i) - lq(i)

				// Penalty:
				-c*gq(i)
				;
	}
	
	/**
	 * Assuming F_i(q_i) = 0.5*q_i^2
	 * => dF_i/dp_i = q_i
	 * @param i
	 * @return
	 */
	private double gradF_q(int i)
	{
		return q(i);
	}
	
	private double gp(int i)
	{
		int dimension = p.getDimension();
		double d = sum(n -> e(i)*e(n)*G(i,n) + f(i)*f(n)*G(i,n) + f(i)*e(n)*B(i,n) - e(i)*f(n)*B(i,n), dimension) - p(i);
		return ZERO_SMALL_G && (d < CON_ERROR_TOLERANCE) ? 0 : d;
	}
	
	private double gq(int i)
	{
		int dimension = q.getDimension();
		double d = sum(n -> f(i)*e(n)*G(i,n) - e(i)*f(n)*G(i,n) - e(i)*e(n)*B(i,n) - f(i)*f(n)*B(i,n), dimension) - q(i);
		return ZERO_SMALL_G && (d < CON_ERROR_TOLERANCE) ? 0 : d;
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
//		int dimension = p.getDimension();
//		double constraintTerms = 
//				 sum(i -> lp(i)*(e(i)*G(i,j) + f(i)*B(i,j)), dimension, j) // i != j
//				+sum(i -> lq(i)*(f(i)*G(i,j) - e(i)*B(i,j)), dimension, j) // i != j
//				+lp(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) // n != j
//					    +2*G(j,j)*e(j))
//				+lq(j)*(sum(n -> -f(n)*G(j,n) - e(n)*B(j,n), dimension, j) // n != j
//					    -2*B(j,j)*e(j));
//		
//		double pt1 = c*sum(i -> i == slackIndex ? 0 : 
//							 gp(i)*(e(i)*G(i,j) + f(i)*B(i,j))
//				   			+gq(i)*(f(i)*G(i,j) - e(i)*B(i,j)),
//				   			dimension, j);
//		double pt2 = +c*gp(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) 
//				+2*e(j)*G(j,j));
//		double pt3 = c*gq(j)*(sum(n -> -f(n)*G(j,n) - e(n)*B(j,n), dimension, j) 
//				-2*e(j)*B(j,j));
//		double penaltyTerms = pt1
//							+
//							(j == slackIndex ? 0 : 
//							(
//								pt2
//								+pt3
//							));
//		
//		return
//			constraintTerms
//			+penaltyTerms
//			;
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
		
		double pt1 = c*sum(i -> i == slackIndex ? 0 : 
						gp(i)*( f(i)*G(i,j) - e(i)*B(i,j))
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
		
//out.print(c*gq(j)+",*(,");
//for(int i = 0; i < dimension; ++i)
//{
//	if(i == 1)
//		out.print(e(i)*G(j,i) - f(i)*B(j,i));
//}
//out.println(",) + "+(-2*f(j)*B(j,j)));
		
		return 
			constraintTerms
			+penaltyTerms
			;
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
	
	protected double lagrangian(RealVector p, RealVector q, RealVector e, RealVector f, RealVector lp, RealVector lq)
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
	
	private double cost(RealVector p, RealVector q)
	{
		int dimension = p.getDimension();
		return 
				sum(i -> generatorMask.getEntry(i) == 1 ? cost_i(p.getEntry(i), q.getEntry(i), p_max.getEntry(i)) : 0, dimension);
//				sum(i -> (1-e(i))*(1-e(i)) + f(i)*f(i), dimension);
	}
	
	private RealVector gp(RealVector p, RealVector q, RealVector e, RealVector f)
	{
		int dimension = p.getDimension();
		return vectorByComponents(i -> gp(p, q, e, f, i), dimension);
	}
	
	private double gp(RealVector p, RealVector q, RealVector e, RealVector f, int i)
	{
		int dimension = p.getDimension();
		double p_i = p.getEntry(i);
		double e_i = e.getEntry(i);
		double f_i = f.getEntry(i);
		double d = sum(n -> e_i*e.getEntry(n)*G(i,n) + f_i*f.getEntry(n)*G(i,n) + f_i*e.getEntry(n)*B(i,n) - e_i*f.getEntry(n)*B(i,n), dimension)
					- p_i;
		return ZERO_SMALL_G && (d < CON_ERROR_TOLERANCE) ? 0 : d;
	}
	
	private RealVector gq(RealVector p, RealVector q, RealVector e, RealVector f)
	{
		int dimension = q.getDimension();
		return vectorByComponents(i -> gq(p, q, e, f, i), dimension);
	}
	
	private double gq(RealVector p, RealVector q, RealVector e, RealVector f, int i)
	{
		int dimension = p.getDimension();
		double q_i = q.getEntry(i);
		double e_i = e.getEntry(i);
		double f_i = f.getEntry(i);
		double d = sum(j -> f_i*e.getEntry(j)*G(i,j) - e_i*f.getEntry(j)*G(i,j) - e_i*e.getEntry(j)*B(i,j) - f_i*f.getEntry(j)*B(i,j), dimension)
					- q_i;
		return ZERO_SMALL_G && (d < CON_ERROR_TOLERANCE) ? 0 : d;
	}

	public static RealVector vectorByComponents(IndexedFunction f, int dimension)
	{
		RealVector grad = new ArrayRealVector(dimension);
		for(int i = 0; i < dimension; ++i)
			grad.setEntry(i, f.value(i));
		
		return grad; 
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
	
	public static void logBusses(AnalysisResults results)
	{
		out.println("index,bus,e,f,p,q");
		Map<String, Integer> numbers = results.getBusNumbers();
		for (String name : numbers.keySet())
		{
			int index = numbers.get(name);
			out.print(index);
			out.print(',');
			out.print(name);
			out.print(',');
			Complex v = results.getBusVoltage(name);
			out.print(format(v.getReal()));
			out.print(',');
			out.print(format(v.getImaginary()));
			out.print(',');
			Complex s = results.getBusPower(name);
			out.print(format(s.getReal()));
			out.print(',');
			out.print(format(s.getImaginary()));
			out.println();
		}
	}
}