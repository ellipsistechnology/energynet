package ellipsis.energy.mas;

import static java.lang.Math.cos;
import static java.lang.Math.sin;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import ellipsis.energy.test.IEEE13BusGrid;

public class PowerFlowAgentProcess implements AgentProcess
{
////Distributed Power Flow ////
	
	private double voltageStepSize = 1e-6;
	private Complex v = Complex.ONE;
	private Complex targetPower;
	private Complex deltaS = Complex.ZERO;
	private double maxPowerError = 1e-3;
	private Agent agent;
	
	private Complex s; // debug only
	
	public PowerFlowAgentProcess()
	{
	}

	/**
	 * Takes an improvement step for the estimation of local voltage.
	 * @return True iff converged.
	 */
	public boolean stepInit()
	{
		if(agent.isSlack())
			return true;

		s = calculatePower();
		deltaS = calculatePowerError(s);

		return deltaS.abs() < maxPowerError;
	}
	
	public void step()
	{
		if(agent.isSlack())
			return;

		// Calculate deltaV from subset of Jacobian:
//		Complex deltaV = calculateVoltageErrors_1();
		
		// Calculate deltaV from full Jacobian with holes:
//		Complex deltaV = calculateVoltageErrors_2();

		// Calculate deltaV from partial, constant sensitivity matrix:
		RealVector deltaV = calculateVoltageErrors_3();
		
		// Calculate new voltage estimate:
		v = addPolar(v, deltaV.mapMultiply(voltageStepSize));
		
//		debug(s, deltaV);
	}

	private Complex addPolar(Complex v1, RealVector v2)
	{
		double abs1 = v1.abs();
		double abs2 = v2.getEntry(1);
		double arg1 = v1.getArgument();
		double arg2 = v2.getEntry(0);
		return ComplexUtils.polar2Complex(
				abs1+abs2,
				arg1+arg2);
	}


	/**
	 * Solves deltaS = J deltaV.
	 * deltaS is a real vector of the power errors and J is the Jacobian matrix.
	 * The Jacobian is the matrix of changes in power with respect to changes in voltage.
	 * @return The local voltage error according to the neighbouring and local power errors.
	 */
	private Complex calculateVoltageErrors_1()
	{
		RealMatrix jacobian = calculateJacobian_1();
		
		SingularValueDecomposition decomp = new SingularValueDecomposition(jacobian);
		DecompositionSolver solver = decomp.getSolver();
		RealVector neighbourDeltaS = neighbourDeltaS();
		RealVector deltaV = solver.solve(neighbourDeltaS);
		return toPolarComplex(deltaV);
	}

	/**
	 * Solves deltaS = J deltaV.
	 * deltaS is a real vector of the power errors and J is the Jacobian matrix.
	 * @param jacobian The jacobian matrix of changes in power with respect to changes in voltage.
	 * @return The local voltage error according to the neighbouring and local power errors.
	 */
	private Complex calculateVoltageErrors_2()
	{
		RealMatrix jacobian = calculateJacobian_2();
		
		SingularValueDecomposition decomp = new SingularValueDecomposition(jacobian);
		DecompositionSolver solver = decomp.getSolver();
		RealVector neighbourDeltaS = neighbourDeltaS();
		RealVector deltaV = solver.solve(neighbourDeltaS);
		
		// Extract local voltage delta:
		int size = deltaV.getDimension()/2;
		double abs = deltaV.getEntry(2*size-1);
		double arg = deltaV.getEntry(size-1);
		return toPolarComplex(abs, arg);
	}
	
	private RealVector calculateVoltageErrors_3()
	{
		RealVector neighbourDeltaS = neighbourDeltaS();
		RealVector deltaV = agent.getSensitivities().operate(neighbourDeltaS);
		return deltaV;
	}

	/**
	 * Combines the neighbours' power errors into a vector of real and imaginary parts
	 * and appends this agent's power error at the end.
	 * @return A real vector of n+1 real values followed by n+1 imaginary values, where n is the number of neighbours.
	 */
	private RealVector neighbourDeltaS()
	{
		int size = agent.getNeighbours().size()+1;
		if(agent.hasSlackNeighbour())
			--size;
		double[] aDeltaS = new double[size*2];
		int j = 0;
		for(Agent agent_j : agent.getNeighbours().keySet())
		{
			if(agent_j.isSlack())
				continue;
			Complex deltaS_j = getDeltaS();
			aDeltaS[j] = deltaS_j.getReal();
			aDeltaS[j+size] = deltaS_j.getImaginary();
			++j;
		}

		aDeltaS[size-1] = deltaS.getReal();
		aDeltaS[2*size-1] = deltaS.getImaginary();
		
		return new ArrayRealVector(aDeltaS);
	}

	/**
	 * Converts a two element real vector into a complex number.
	 * @param v A vector with two entries: The complex argument followed by the absolute value.
	 * @return A complex representation of the given vector.
	 */
	private Complex toPolarComplex(RealVector v)
	{
		double abs = v.getEntry(1);
		double arg = v.getEntry(0);
		return toPolarComplex(abs, arg);
	}

	public Complex toPolarComplex(double abs, double arg)
	{
		if(abs < 0)
		{
			return ComplexUtils.polar2Complex(-abs, arg+Math.PI);
		}
		else
		{
			return ComplexUtils.polar2Complex(abs, arg);
		}
	}

	private RealMatrix calculateJacobian_1()
	{
		int size = agent.getNeighbours().size()+1; // Include self.
		if(agent.hasSlackNeighbour())
			--size;
		double[][] jacobian = new double[2*size][2];
		double v_j_abs = v.abs();
		double v_j_arg = v.getArgument();
		
		// Neighbours:
		{
			int i = 0;
			for (Agent agent_i : agent.getNeighbours().keySet())
			{
				if(agent_i.isSlack())
					continue;
				Complex y_ij = agent.getNeighbours().get(agent_i);
				double y_ij_abs = y_ij.abs();
				double y_ij_arg = y_ij.getArgument();
				Complex v_i = ((PowerFlowAgentProcess)agent_i.getAgentProcess()).getV();
				double v_i_abs = v_i.abs();
				double v_i_arg = v_i.getArgument();
				jacobian[i][0] = dP_idd_j(y_ij_abs, y_ij_arg, v_i_abs, v_i_arg, v_j_abs, v_j_arg); // dP/dd (top-left):
				jacobian[i][1] = dP_idv_j(y_ij_abs, y_ij_arg, v_i_abs, v_i_arg, v_j_arg); // dP/dv (top-right):
				jacobian[i + size][0] = dQ_idd_j(y_ij_abs, y_ij_arg, v_i_abs, v_i_arg, v_j_abs, v_j_arg); // dQ/dd (bottom-left):
				jacobian[i + size][1] = dQ_idv_j(y_ij_abs, y_ij_arg, v_i_abs, v_i_arg, v_j_abs, v_j_arg); // dQ/dv (bottom-right):
						
				++i;
			}
		}
		
		// Self:
		double y_ii_abs = agent.getSelfAdmittance().abs();
		double y_ii_arg = agent.getSelfAdmittance().getArgument();
		jacobian[size-1][0] = dP_idd_i(v_j_abs, v_j_arg); // dP/dd (top-left):
		jacobian[size-1][1] = dP_idv_i(y_ii_abs, y_ii_arg, v_j_abs, v_j_arg); // dP/dv (top-right):
		jacobian[2*size-1][0] = dQ_idd_i(v_j_abs, v_j_arg); // dQ/dd (bottom-left):
		jacobian[2*size-1][1] = dQ_idv_i(y_ii_abs, y_ii_arg, v_j_abs, v_j_arg); // dQ/dv (bottom-right):
		
		return new Array2DRowRealMatrix(jacobian);
	}

	private RealMatrix calculateJacobian_2()
	{
		int size = agent.getNeighbours().size()+1; // Include self.
		if(agent.hasSlackNeighbour())
			--size;
		double[][] jacobian = new double[2*size][2*size];
		
		// Neighbours i:
		{
			int i = 0;
			for (Agent agent_i : agent.getNeighbours().keySet())
			{
				if(agent_i.isSlack())
					continue;

				// Get agent i's voltage:
				Complex v_i = ((PowerFlowAgentProcess)agent_i.getAgentProcess()).getV();
				double v_i_abs = v_i.abs();
				double v_i_arg = v_i.getArgument();
				
				// Neighbours j:
				{
					int j = 0;
					for (Agent agent_j : agent.getNeighbours().keySet())
					{
						fillInJacobian(size, jacobian, i, agent_i, v_i_abs, v_i_arg, j, agent_j);
						++j;
					}
				}
				
				// From neighbour i to self:
				fillInJacobian(size, jacobian, i, agent_i, v_i_abs, v_i_arg, size-1, agent);
						
				++i;
			}
		}
		
		// Neighbours j:
		{
			int j = 0;
			double v_i_abs = v.abs();
			double v_i_arg = v.getArgument();
			for (Agent agent_j : agent.getNeighbours().keySet())
			{
				// From self to neighbour j:
				fillInJacobian(size, jacobian, size-1, agent, v_i_abs, v_i_arg, j, agent_j);
				++j;
			}
			
			// From self to self:
			fillInJacobian(size, jacobian, size-1, agent, v_i_abs, v_i_arg, size-1, agent);
		}
		
		return new Array2DRowRealMatrix(jacobian);
	}


	public void fillInJacobian(int size, double[][] jacobian, int i,
			Agent agent_i, double v_i_abs, double v_i_arg, int j, Agent agent_j)
	{
		// Get agent j's voltage:
		double v_j_abs = v.abs();
		double v_j_arg = v.getArgument();

		// Get admittance between j and i:
		// FIXME This is cheating at the moment:
		Complex y_ij;
		if(i == j)
			y_ij = agent_i.getSelfAdmittance();
		else
			y_ij = agent_j.getNeighbours().get(agent_i);
		if(y_ij == null) // For radial systems this will always be the case - so maybe not cheating after all.
			y_ij = Complex.ZERO;
		double y_ij_abs = y_ij.abs();
		double y_ij_arg = y_ij.getArgument();
		
		// Calculate off-diagonel Jacobian elements:
		if(i != j)
		{
			jacobian[i][j] = dP_idd_j(y_ij_abs, y_ij_arg, v_i_abs, v_i_arg, v_j_abs, v_j_arg); // dP/dd (top-left):
			jacobian[i][j+size] = dP_idv_j(y_ij_abs, y_ij_arg, v_i_abs, v_i_arg, v_j_arg); // dP/dv (top-right):
			jacobian[i+size][j] = dQ_idd_j(y_ij_abs, y_ij_arg, v_i_abs, v_i_arg, v_j_abs, v_j_arg); // dQ/dd (bottom-left):
			jacobian[i+size][j+size] = dQ_idv_j(y_ij_abs, y_ij_arg, v_i_abs, v_i_arg, v_j_abs, v_j_arg); // dQ/dv (bottom-right):
		}
		else
		{
			jacobian[i][j] = dP_idd_i(v_j_abs, v_j_arg); // dP/dd (top-left):
			jacobian[i][j+size] = dP_idv_i(y_ij_abs, y_ij_arg, v_j_abs, v_j_arg); // dP/dv (top-right):
			jacobian[i+size][j] = dQ_idd_i(v_j_abs, v_j_arg); // dQ/dd (bottom-left):
			jacobian[i+size][j+size] = dQ_idv_i(y_ij_abs, y_ij_arg, v_j_abs, v_j_arg); // dQ/dv (bottom-right):
		}
	}

	private double dP_idd_j(double y_ij_abs, double Y_ij_arg, double v_i_abs, double v_i_arg, double v_j_abs, double v_j_arg) 
	{
		return -v_i_abs*v_j_abs*y_ij_abs*sin(Y_ij_arg - v_i_arg + v_j_arg);
	}
	
	private double dP_idv_j(double Y_ij_abs, double Y_ij_arg, double v_i_abs, double v_i_arg, double v_j_arg) 
	{
		return v_i_abs*Y_ij_abs*cos(Y_ij_arg - v_i_arg + v_j_arg);
	}

	private double dQ_idv_j(double Y_ij_abs, double Y_ij_arg, double v_i_abs, double v_i_arg, double v_j_abs, double v_j_arg) 
	{
		return -v_i_abs*Y_ij_abs*sin(Y_ij_arg - v_i_arg + v_j_arg);
	}

	private double dQ_idd_j(double Y_ij_abs, double Y_ij_arg, double v_i_abs, double v_i_arg, double v_j_abs, double v_j_arg) 
	{
		return -v_i_abs*v_j_abs*Y_ij_abs*cos(Y_ij_arg - v_i_arg + v_j_arg);
	}

	private double dP_idd_i(double v_i_abs, double v_i_arg) 
	{
		double d = 0;
		for (Agent agent_j : agent.getNeighbours().keySet())
		{
			Complex y_ij = agent.getNeighbours().get(agent_j);
			double y_ij_abs = y_ij.abs();
			double y_ij_arg = y_ij.getArgument();
			Complex v_j = ((PowerFlowAgentProcess)agent_j.getAgentProcess()).getV();
			double v_j_abs = v_j.abs();
			double v_j_arg = v_j.getArgument();
			d += v_i_abs*v_j_abs*y_ij_abs*sin(y_ij_arg - v_i_arg + v_j_arg);
		}
		return d;
	}

	private double dP_idv_i(double y_ii_abs, double y_ii_arg, double v_i_abs, double v_i_arg) 
	{
		double d = 0;
		for (Agent agent_j : agent.getNeighbours().keySet())
		{
			Complex y_ij = agent.getNeighbours().get(agent_j);
			double y_ij_abs = y_ij.abs();
			double y_ij_arg = y_ij.getArgument();
			Complex v_j = ((PowerFlowAgentProcess)agent_j.getAgentProcess()).getV();
			double v_j_abs = v_j.abs();
			double v_j_arg = v_j.getArgument();
			d += v_j_abs*y_ij_abs*cos(y_ij_arg - v_i_arg + v_j_arg);
		}
		return 2*v_i_abs*y_ii_abs*cos(y_ii_arg) + d;
	}

	private double dQ_idd_i(double v_i_abs, double v_i_arg) 
	{
		double d = 0;
		for (Agent agent_j : agent.getNeighbours().keySet())
		{
			Complex y_ij = agent.getNeighbours().get(agent_j);
			double y_ij_abs = y_ij.abs();
			double y_ij_arg = y_ij.getArgument();
			Complex v_j = ((PowerFlowAgentProcess)agent_j.getAgentProcess()).getV();
			double v_j_abs = v_j.abs();
			double v_j_arg = v_j.getArgument();
			return v_i_abs*v_j_abs*y_ij_abs*cos(y_ij_arg - v_i_arg + v_j_arg);
		}
		
		return d;
	}

	private double dQ_idv_i(double y_ii_abs, double y_ii_arg, double v_i_abs, final double v_i_arg) 
	{
		double d = 0;
		for (Agent agent_j : agent.getNeighbours().keySet())
		{
			Complex y_ij = agent.getNeighbours().get(agent_j);
			double y_ij_abs = y_ij.abs();
			double y_ij_arg = y_ij.getArgument();
			Complex v_j = ((PowerFlowAgentProcess)agent_j.getAgentProcess()).getV();
			double v_j_abs = v_j.abs();
			double v_j_arg = v_j.getArgument();
			d += v_j_abs*y_ij_abs*sin(y_ij_arg - v_i_arg + v_j_arg);
		}
		
		return -2*v_i_abs*y_ii_abs*sin(y_ii_arg) - d;
	}

	private Complex calculatePowerError(Complex s)
	{
		return targetPower.subtract(s);
	}

	private Complex calculatePower()
	{
//		double p = 0;
//		double q = 0;
		Complex s = Complex.ZERO;
//		double v_i_abs = v.abs();
//		double v_i_arg = v.getArgument();
		
		// Neighbours:
		for (Agent agent_j : agent.getNeighbours().keySet())
		{
			Complex y_ij = agent.getNeighbours().get(agent_j);
			Complex v_j = ((PowerFlowAgentProcess)agent_j.getAgentProcess()).getV();
//			p += v_i_abs*v_j.abs()*y_ij.abs()*cos(y_ij.getArgument()-v_i_arg+v_j.getArgument());
//			q += -v_i_abs*v_j.abs()*y_ij.abs()*sin(y_ij.getArgument()-v_i_arg+v_j.getArgument());
			s = s.add(y_ij.multiply(v_j));
		}
		
		// Self:
//		p += v_i_abs*v_i_abs*selfAdmittance.abs()*cos(selfAdmittance.getArgument());
//		q += -v_i_abs*v_i_abs*selfAdmittance.abs()*sin(selfAdmittance.getArgument());
		s = s.add(agent.getSelfAdmittance().multiply(v));
		
		//
		s = s.multiply(v.conjugate()).conjugate(); // (v_i* sum{y_ij v_j})*
		
//System.out.println(bus.getName()+" ["+p+" + j"+q+"] = "+s);
		
		return s;
	}
	
	
	//// Accessors ////
	
	public void setAgent(Agent agent)
	{
		this.agent = agent;
		targetPower = agent.getBus().getNetPower(false).divide(IEEE13BusGrid.BASE_POWER); // FIXME
	}

	public Complex getV()
	{
		return v;
	}
	
	public Complex getDeltaS()
	{
		return deltaS;
	}
	
	
	//// Testing ////

	public void debug(Complex s, RealVector deltaV)
	{
//		if(bus.getName().equals("684"))
//			System.out.println(
//							s.getReal()+","
//							+s.getImaginary()+","
//							+targetPower.getReal()+","
//							+targetPower.getImaginary()+","
//							+deltaS.getReal()+","
//							+deltaS.getImaginary()+","
//							+v.abs()+","
//							+v.getArgument()+","
//							+deltaV.abs()+","
//							+deltaV.getArgument()
//					);
		if(agent.getBus().getName().equals("671"))
		{
			System.out.print(
							 v.abs()+","
							+v.getArgument()+","
							+deltaV.getEntry(1)+","
							+deltaV.getEntry(0)+","
							+s.getReal()+","
							+s.getImaginary()+","
							+targetPower.getReal()+","
							+targetPower.getImaginary()+","
							+deltaS.getReal()+","
							+deltaS.getImaginary()+","
					);
			for (Agent agent_j : agent.getNeighbours().keySet())
			{
				Complex deltaS_j = ((PowerFlowAgentProcess)agent_j.getAgentProcess()).getDeltaS();
				System.out.print(deltaS_j.getReal()+",");
			}
			for (Agent agent_j : agent.getNeighbours().keySet())
			{
				Complex deltaS_j = ((PowerFlowAgentProcess)agent_j.getAgentProcess()).getDeltaS();
				System.out.print(deltaS_j.getImaginary()+",");
			}
			System.out.println();
		}
	}
}
