package ellipsis.energy.scenario;

import ellipsis.energy.calculation.AnalysisResults;

public interface Regulator
{
    void regulate(AnalysisResults results);
}
