package ellipsis.energy.calculation;

import org.apache.commons.math3.complex.Complex;

public class GaussSeidelSolver implements PowerFlowSolver
{
    private boolean[] isSlackSource;
    private Complex[] power;
    private double error;
    private AdmittanceMatrix admittances;
    
    public GaussSeidelSolver(Complex[] power, boolean[] isSlackSource, AdmittanceMatrix admittances)
    {
        this.power = power;
        this.isSlackSource = isSlackSource;
        this.admittances = admittances;
    }

    @Override
    public Complex[] solve(Complex[] Vk)
    {
        error = 0;
        int nBusses = Vk.length;
        
        for (int i = 0; i < nBusses; ++i)
        {
            if(isSlackSource[i])
                continue;
            
            // sum(Y_ij*Vj(k), i != j):
            Complex vSum = Complex.ZERO;
            for (int j = 0; j < nBusses; ++j)
            {
                if(i != j)
                {
                    vSum = vSum.add( admittances.get(i, j).multiply(Vk[j]) );
                }
            }
            
            // V_i(k+1) = {S_i/V_i*(k) + sum(y_ij*Vj(k), i != j)}/sum(y_ij)
            Complex Vkplus1 = power[i].conjugate().divide(Vk[i].conjugate()).subtract(vSum).divide(admittances.get(i, i));
/*System.out.println(String.format(
      "%d:%d: %s = {%s/%s + %s}/%s",
      k,
      i,
      Vkplus1.toString(),
      power[i].conjugate().toString(),
      Vk[i].conjugate().toString(),
      vSum.negate(),
      admittances.getEntry(i, i)));*/
            
            // Get the maximum error:
            double newError = Math.abs(Vk[i].subtract(Vkplus1).abs());
            if(newError > error)
                error = newError;
//if(i == 0) System.out.println("iteration "+k+", bus at index 0 has error "+error);
            
            Vk[i] = Vkplus1;
        }
        
        return Vk;
    }

    @Override
    public double getError()
    {
        return error;
    }
}
