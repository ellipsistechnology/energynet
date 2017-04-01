package ellipsis.util;

import la.matrix.DenseMatrix;
import la.matrix.Matrix;
import ml.utils.InPlaceOperator;
import ml.utils.Matlab;

/**
 * Solves a constrained problem of the form:
 * Minimise  f(x),
 * Such that F_i(x) < 0, for all i,
 * And       Ax = b.
 * 
 * @author bmillar
 * @author Mingjie Qian
 */
public class PrimalDualInteriorPointSolver
{
	private boolean verbose = false;
	
	public static interface Constraints
	{
		Matrix derivative(Matrix x);
		Matrix value(Matrix x);
		Matrix hessian(int i, Matrix x);
		int constrantCount();
	}
	
	public static interface Objective
	{
		double value(Matrix x);
		Matrix hessian(Matrix x);
		Matrix gradient(Matrix x);
	}

	private Objective f;
	private Constraints F;
	private Matrix A;
	private Matrix b;
	private int maxIterations = 1000;
//	private boolean equalityConstraintIsZero;
    
    public PrimalDualInteriorPointSolver()
    {
    }

    /**
     * 
     * @param f The objective function.
     * @param F The inequality constraints.
     * @param A The linear equality matrix: Ax = b 
     * @param b The linear equality vector: Ax = b
     */
    public PrimalDualInteriorPointSolver(Objective f, Constraints F, Matrix A, Matrix b)
	{
		this.f = f;
		this.F = F;
		this.A = A;
		this.b = b;
		
//		equalityConstraintIsZero = true;
//		for(int i = 0; i < b.getRowDimension(); ++i)
//		{
//			if(b.getEntry(i, 0) != 0)
//			{
//				equalityConstraintIsZero = false;
//				break;
//			}
//		}
	}

public static int convergedIterations = 0;
public static int nonconvergedIterations = 0;
	public boolean solve(Matrix x0)
	{
		// Init:
		Matrix x = x0.copy();
		init(x);
		double fval = f.value(x);
		Matrix DF_x = F.derivative(x);
		Matrix F_x = F.value(x);
		Matrix G_f_x = f.gradient(x);
		Matrix H_x = hessian(x, l);
		
		// To be filled in by step():
		Matrix l = new DenseMatrix(F_x.getRowDimension(), 1);
		Matrix v = new DenseMatrix(A.getRowDimension(), 1);
		
		// Solve:
		boolean flags[] = null;
		int k = 0;
		while (k < maxIterations)
		{
//verboseln(k+":");
//verboseln("x="+x);
//verboseln("G_f_x="+G_f_x);
//verboseln("H_x="+H_x);
			flags = step(A, b, H_x, F_x, DF_x, G_f_x, fval, x, l, v);
			if (flags[0])
				break;
			// Compute the objective function value, if flags[1] is true
			// gradient will also be computed.
			fval = f.value(x);
			F_x = F.value(x);
//			if (flags[1])
			{
				// Compute the gradient
				G_f_x = f.gradient(x);
				H_x = hessian(x, l);
				DF_x = F.derivative(x);
			}
			
			++k;
		}
if(flags[0])
	++convergedIterations;
else
	++nonconvergedIterations;

		state = 0;
		
		return flags[0];
	}

	public void init(Matrix x)
	{
        la.matrix.Matrix F_x = F.value(x);
		this.n = A.getColumnDimension();
		this.p = A.getRowDimension();
		this.m = F_x.getRowDimension();
        this.x = x.copy();
        this.l = Matlab.rdivide(Matlab.ones(m, 1), m);
        this.v = Matlab.zeros(p, 1);
        this.l_s = l.copy();
        this.v_s = v.copy();
        this.eta_t = -Matlab.innerProduct(F_x, l);
	}

	private la.matrix.Matrix hessian(Matrix x, Matrix l)
	{
		la.matrix.Matrix Hf = f.hessian(x);
		
		int dim = x.getRowDimension();
		Matrix sum = new DenseMatrix(dim, dim);
		for(int i = 0; i < F.constrantCount(); ++i)
		{
			double l_i = l.getEntry(i, 0);
			Matrix summand = F.hessian(i, x).times(l_i);
			InPlaceOperator.plusAssign(sum, summand);
		}
		return Hf.plus(sum);
	}
	
	
	//// Internal Variables and Process ////
	
    private Matrix x = null;
    private Matrix l = null;
    private Matrix v = null;
    private boolean gradientRequired = false;
    private boolean converge = false;
    private int state = 0;
    private double t = 1.0D;
    private int n = 0;
    private int p = 0;
    private int m = 0;
    private double mu = 10;//1.8D;
    private double epsilon = 1e-5;//1E-010D;
    private double epsilon_feas = 1e-5;//1E-010D;
    private double alpha = 0.01;//0.10000000000000001D;
    private double beta = 0.5;//0.97999999999999998D;
    private double eta_t = 1.0D;
    private double s = 1.0D;
    private double residual = 0.0D;
    private Matrix r_prim = null;
    private Matrix r_dual = null;
    private Matrix r_cent = null;
    private Matrix Matrix = null;
    private Matrix Vector = null;
    private Matrix z_pd = null;
    private Matrix x_nt = null;
    private Matrix l_nt = null;
    private Matrix v_nt = null;
    private Matrix l_s = null;
    private Matrix v_s = null;

	public Matrix getOptimalLambda()
    {
        return l;
    }

    public Matrix getOptimalNu()
    {
        return v;
    }

    public boolean[] step(Matrix A, Matrix b, Matrix H_x, Matrix F_x, Matrix DF_x, Matrix G_f_x, double fval, Matrix x)
    {
        return step(A, b, H_x, F_x, DF_x, G_f_x, fval, x, null, null);
    }

    public boolean[] step(Matrix A, Matrix b, Matrix H_x, Matrix F_x, Matrix DF_x, Matrix G_f_x, double fval, Matrix x_s, Matrix l_opt, Matrix v_opt)
    {
        if(state == 5)
        {
verboseln("state 5");
            state = 0;
        }
        if(state == 0)
        {
verboseln("state 0");
            state = 1;
        }
        if(state == 1)
        {
verboseln("state 1");
            double residual_prim = 0.0D;
            double residual_dual = 0.0D;
            t = (mu * (double)m) / eta_t;
            r_prim = A.mtimes(x).minus(b);
            r_dual = G_f_x.plus(DF_x.transpose().mtimes(l)).plus(A.transpose().mtimes(v));
            r_cent = Matlab.uminus(Matlab.times(l, F_x)).minus(Matlab.rdivide(Matlab.ones(m, 1), t));
            Matrix = Matlab.vertcat(new Matrix[] {
                Matlab.horzcat(new Matrix[] {
                    H_x, DF_x.transpose(), A.transpose()
                }), Matlab.horzcat(new Matrix[] {
                    Matlab.uminus(Matlab.mtimes(Matlab.diag(l), DF_x)), Matlab.uminus(Matlab.diag(F_x)), Matlab.zeros(m, p)
                }), Matlab.horzcat(new Matrix[] {
                    A, Matlab.zeros(p, m), Matlab.zeros(p, p)
                })
            });
            Vector = Matlab.uminus(Matlab.vertcat(new Matrix[] {
                r_dual, r_cent, r_prim
            }));
            residual = Matlab.norm(Vector);
            residual_prim = Matlab.norm(r_prim);
            residual_dual = Matlab.norm(r_dual);
            eta_t = -Matlab.innerProduct(F_x, l);
verboseln("resuduals: "+residual_prim+", "+residual_dual+", "+eta_t);
            if(residual_prim <= epsilon_feas && residual_dual <= epsilon_feas && eta_t <= epsilon)
            {
//                Printer.fprintf("Terminate successfully.\n\n", new Object[0]);
                if(l_opt != null)
                    Matlab.setMatrix(l_opt, l);
                if(v_opt != null)
                    Matlab.setMatrix(v_opt, v);
                converge = true;
                gradientRequired = false;
                state = 5;
//                System.out.printf("Primal-dual interior-point algorithm converges.\n", new Object[0]);
                return (new boolean[] {
                    converge, gradientRequired
                });
            }
            z_pd = Matlab.mldivide(Matrix, Vector);//solveLinearMatrix(Matrix, Vector);
            x_nt = Matlab.getRows(z_pd, 0, n - 1);
            l_nt = Matlab.getRows(z_pd, n, (n + m) - 1);
            v_nt = Matlab.getRows(z_pd, n + m, (n + m + p) - 1);
            s = 1.0D;
            do
            {
                InPlaceOperator.affine(l_s, s, l_nt, '+', l);
                if(Matlab.sumAll(Matlab.lt(l_s, 0.0D)) > 0.0D)
                {
                    s = beta * s;
                } else
                {
                    state = 2;
                    InPlaceOperator.affine(x_s, s, x_nt, '+', x);
                    converge = false;
                    gradientRequired = false;
                    return (new boolean[] {
                        converge, gradientRequired
                    });
                }
            } while(true);
        }
        if(state == 2)
        {
verboseln("state 2");
            if(Matlab.sumAll(Matlab.gt(F_x, 0.0D)) > 0.0D)
            {
                s = beta * s;
                InPlaceOperator.affine(x_s, s, x_nt, '+', x);
                converge = false;
                gradientRequired = false;
                return (new boolean[] {
                    converge, gradientRequired
                });
            } else
            {
                state = 3;
                converge = false;
                gradientRequired = true;
                return (new boolean[] {
                    converge, gradientRequired
                });
            }
        }
        if(state == 3)
        {
            Matrix r_prim_s = null;
            Matrix r_dual_s = null;
            Matrix r_cent_s = null;
            double residual_s = 0.0D;
            InPlaceOperator.affine(l_s, s, l_nt, '+', l);
            InPlaceOperator.affine(v_s, s, v_nt, '+', v);
            r_prim_s = A.mtimes(x_s).minus(b);
            r_dual_s = G_f_x.plus(DF_x.transpose().mtimes(l_s)).plus(A.transpose().mtimes(v_s));
            r_cent_s = Matlab.uminus(Matlab.times(l_s, F_x)).minus(Matlab.rdivide(Matlab.ones(m, 1), t));
            residual_s = Matlab.norm(Matlab.vertcat(new Matrix[] {
                r_dual_s, r_cent_s, r_prim_s
            }));

verboseLog(x_s, r_prim_s, r_dual_s, r_cent_s, residual_s);
            if(residual_s <= (1.0D - alpha * s) * residual)
            {
                InPlaceOperator.assign(l, l_s);
                InPlaceOperator.assign(v, v_s);
                InPlaceOperator.assign(x, x_s);
                state = 4;
            } else
            {
                s = beta * s;
                converge = false;
                gradientRequired = true;
                InPlaceOperator.affine(x_s, s, x_nt, '+', x);
                return (new boolean[] {
                    converge, gradientRequired
                });
            }
        }
        if(state == 4)
        {
verboseln("state 4");
            state = 1;
        }
        converge = false;
        gradientRequired = true;
        return (new boolean[] {
            converge, gradientRequired
        });
    }

	private void verboseLog(Matrix x_s, Matrix r_prim_s, Matrix r_dual_s,
			Matrix r_cent_s, double residual_s)
	{
		if(!verbose)
			return;
		verboseln("state 3, x=["+x_s.getEntry(0, 0)+","+x_s.getEntry(1, 0)+"], \tr(x,l,v)=["
				+ r_dual_s.getEntry(0, 0)
				+ ","
				+ r_dual_s.getEntry(1, 0)
				+ ","
				+ r_cent_s.getEntry(0, 0)
				+ ","
//				+ r_cent_s.getEntry(1, 0)
//				+ ","
				+ r_prim_s.getEntry(0, 0)
				+ "] \t"+residual_s+" <= (1.0D - "+alpha+" * "+s+") * "+residual+"="+((1.0D - alpha * s) * residual));
	}

//	private la.matrix.Matrix solveLinearMatrix(Matrix m, Matrix v)
//	{
//		int ARowDimension = 0;
//		if(equalityConstraintIsZero)
//		{
//			ARowDimension = A.getRowDimension();
//			m = m.getSubMatrix(0, m.getRowDimension()-ARowDimension-1, 0, m.getColumnDimension()-ARowDimension-1); // startRow, endRow, startCol, endCol
//			// Note that both use A.rowDimension because one is transposed.
//			v = v.getSubMatrix(0, v.getRowDimension()-ARowDimension-1, 0, 0);
//		}
//		Matrix solution;
//		try
//		{
//			RealMatrix rm = new Array2DRowRealMatrix(m.getData());
//			RealVector rv = new ArrayRealVector(v.getColumnVector(0).getData());
//			RealVector res = new QRDecomposition(rm).getSolver().solve(rv);
//			solution = toMatrix(res);
//		}
//		catch(SingularMatrixException e)
//		{
//			solution = Matlab.mldivide(m, v);
//		}
//		
//		if(equalityConstraintIsZero)
//		{
//			int dim = solution.getRowDimension();
//			Matrix fullSolution = new DenseMatrix(dim+ARowDimension, 1);
//			for(int i = 0; i < dim; ++i)
//			{
//				fullSolution.setEntry(i, 0, solution.getEntry(i, 0));
//			}
//			solution = fullSolution;
//		}
//		
//		return solution;
//	}

//	private Matrix toMatrix(RealVector v)
//	{
//		int dimension = v.getDimension();
//		double[][] m = new double[dimension][1];
//		for(int i = 0; i < dimension; ++i)
//		{
//			m[i][0] = v.getEntry(i);
//		}
//		return new DenseMatrix(m);
//	}
    
    
    //// Debug ////

    protected void verboseln(String string)
	{
		if(verbose) System.out.println(string);
	}

	protected void verbose(String string)
	{
		if(verbose) System.out.print(string);
	}
	
	
	//// Accessors ////

	public int getMaxIterations()
	{
		return maxIterations;
	}

	public void setMaxIterations(int maxIterations)
	{
		this.maxIterations = maxIterations;
	}

	public Matrix getX()
	{
		return x;
	}

	public Objective getf()
	{
		return f;
	}

	public void setf(Objective f)
	{
		this.f = f;
	}

	public Constraints getF()
	{
		return F;
	}

	public void setF(Constraints F)
	{
		this.F = F;
	}

	public Matrix getA()
	{
		return A;
	}

	public void setA(Matrix A)
	{
		this.A = A;
	}

	public Matrix getb()
	{
		return b;
	}

	public void setb(Matrix b)
	{
		this.b = b;
	}

	public double getEpsilon()
	{
		return epsilon;
	}

	public void setEpsilon(double epsilon)
	{
		this.epsilon = epsilon;
	}

	public double getEpsilon_feas()
	{
		return epsilon_feas;
	}

	public void setEpsilon_feas(double epsilon_feas)
	{
		this.epsilon_feas = epsilon_feas;
	}

	public double getAlpha()
	{
		return alpha;
	}

	public void setAlpha(double alpha)
	{
		this.alpha = alpha;
	}

	public double getBeta()
	{
		return beta;
	}

	public void setBeta(double beta)
	{
		this.beta = beta;
	}
}
