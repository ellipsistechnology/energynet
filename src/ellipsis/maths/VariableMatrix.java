package ellipsis.maths;

public class VariableMatrix
{
	private Variable[][] vars;
	
	public VariableMatrix(int rows, int cols)
	{
		if(rows <= 0 || cols <= 0)
			throw new RuntimeException("Row and column dimensions must be positive.");
		vars = new Variable[rows][cols];
	}
	
	public VariableMatrix(Variable[][] vars)
	{
		this.vars = vars;
	}
	
	public Variable getElementAt(int row, int col)
	{
		Variable v = vars[row][col];
		if(v == null)
			return Variable.ZERO;
		return v;
	}
	
	public void setElementAt(int row, int col, Variable v)
	{
		vars[row][col] = v;
	}
	
	public int getRowDimension()
	{
		return vars.length;
	}
	
	public int getColumnDimension()
	{
		return vars[0].length;
	}
	
	public VariableMatrix multiply(VariableMatrix m)
	{
		int rowDimension = getRowDimension();
		int columnDimension = getColumnDimension();
		if(m.getRowDimension() != columnDimension)
			throw new RuntimeException("Given matrix row dimension did not match this matrix's column dimension.");
		
		Variable[][] product = new Variable[m.getRowDimension()][m.getColumnDimension()];
		
		for(int i = 0; i < rowDimension; ++i)
		{
			for(int j = 0; j < columnDimension; ++j)
			{
				product[i][j] = multiply(vars[i], m, j);
			}
		}
		
		return new VariableMatrix(product);
	}

	private Variable multiply(Variable[] row, VariableMatrix m, int colIndex)
	{
		Variable v = Variable.ZERO;
		for(int i = 0; i < row.length; ++i)
		{
			Variable m_ij = m.getElementAt(i, colIndex);
			v = v.add(row[i].multiply(m_ij));
		}
		return v;
	}
	
	@Override
	public String toString()
	{
		StringBuffer sb = new StringBuffer();
		sb.append('[');
		for(int i = 0; i < getRowDimension(); ++i)
		{
			sb.append('[');
			for(int j = 0; j < getColumnDimension(); ++j)
			{
				sb.append(vars[i][j]);
				sb.append(',');
			}
			sb.append("],\n");
		}
		sb.append(']');
		
		return sb.toString();
	}
	
	public String toValueString()
	{
		StringBuffer sb = new StringBuffer();
		sb.append('[');
		for(int i = 0; i < getRowDimension(); ++i)
		{
			sb.append('[');
			for(int j = 0; j < getColumnDimension(); ++j)
			{
				sb.append(vars[i][j].getValue());
				sb.append(',');
			}
			sb.append("],\n");
		}
		sb.append(']');
		
		return sb.toString();
	}
	
	public VariableMatrix transpose()
	{
		Variable[][] data = new Variable[getRowDimension()][getColumnDimension()];
		for(int i = 0; i < getRowDimension(); ++i)
		{
			for(int j = 0; j < getColumnDimension(); ++j)
			{
				data[j][i] = this.vars[i][j];
			}
		}
		return new VariableMatrix(data);
	}
	
	
	//// Test ////
	
	// FIXME Warning: This is giving incorrect values - anything could be wrong!
	public static void main(String[] args)
	{
		Variable[][] vars = new Variable[][]{
				{new Variable("A11"),new Variable("A12")},
				{new Variable("A21"),new Variable("A22")}
		};
		
		VariableMatrix A = new VariableMatrix(vars);
		A.vars[0][0].setValue(0.49);
		A.vars[0][1].setValue(0.49);
		A.vars[1][0].setValue(0.49);
		A.vars[1][1].setValue(0.49);
		
		System.out.println("A=");
		System.out.println(A);
		System.out.println(A.toValueString());
		System.out.println();
		
		System.out.println("AA=");
		VariableMatrix AA = A.multiply(A);
		System.out.println(AA);
		System.out.println(AA.toValueString());
		System.out.println();
		
		System.out.println("AAA=");
		VariableMatrix AAA = AA.multiply(A);
		System.out.println(AAA);
		System.out.println(AAA.toValueString());
		System.out.println();
		
		System.out.println("AAAA=");
		VariableMatrix AAAA = AAA.multiply(A);
		System.out.println(AAAA);
		System.out.println(AAAA.toValueString());
		System.out.println();
		
		System.out.println("AAAAA=");
		VariableMatrix AAAAA = AAAA.multiply(A);
		System.out.println(AAAAA);
		System.out.println(AAAAA.toValueString());
		System.out.println();
		
		vars = new Variable[][]{
				{new Variable("B11"),new Variable("B12")},
				{new Variable("B21"),new Variable("B22")}
		};
		
		VariableMatrix B = new VariableMatrix(vars);
		System.out.println(A.multiply(B));
		
		System.out.println("AB^T=\n"+A.multiply(B.transpose()));
	}
}
