package ellipsis.energy.calculation;

import org.apache.commons.math3.complex.Complex;

public interface PowerFlowSolver
{
    Complex[] solve(Complex[] Vk);
    double getError();
}
