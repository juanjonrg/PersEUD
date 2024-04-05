package org.uma.jmetal.problem;

import java.io.Serializable;
import java.util.List;

import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.solution.doublesolution.DoubleSolution;

/**
 * Interface representing a multi-objective optimization problem
 *
 * @author Antonio J. Nebro <antonio@lcc.uma.es>
 *
 * @param <S> Encoding
 */
public interface Problem<S> extends Serializable {
  /* Getters */
  int getNumberOfVariables() ;
  int getNumberOfObjectives() ;
  int getNumberOfConstraints() ;
  String getName() ;

  /* Methods */
  void evaluate(S solution) ;
  void evaluate(List<DoubleSolution> solutionList) ;
  S createSolution() ;
}
