package org.mine.exec;

import java.io.FileNotFoundException;
import java.util.List;

import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.algorithm.multiobjective.moead.AbstractMOEAD;
import org.uma.jmetal.algorithm.multiobjective.moead.MOEADBuilder;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.impl.crossover.DifferentialEvolutionCrossover;
import org.uma.jmetal.operator.impl.crossover.SBXCrossover;
import org.uma.jmetal.operator.impl.mutation.PolynomialMutation;
import org.uma.jmetal.problem.DoubleProblem;
import org.uma.jmetal.problem.multiobjective.cdtlz.C1_DTLZ1;
import org.uma.jmetal.problem.multiobjective.cdtlz.C1_DTLZ3;
import org.uma.jmetal.problem.multiobjective.cdtlz.C2_DTLZ2;
import org.uma.jmetal.problem.multiobjective.cdtlz.C3_DTLZ1;
import org.uma.jmetal.problem.multiobjective.cdtlz.C3_DTLZ4;
import org.uma.jmetal.problem.multiobjective.cdtlz.ConvexC2_DTLZ2;
import org.uma.jmetal.problem.multiobjective.dtlz.Convex_DTLZ2;
import org.uma.jmetal.problem.multiobjective.dtlz.DTLZ1;
import org.uma.jmetal.problem.multiobjective.dtlz.DTLZ2;
import org.uma.jmetal.problem.multiobjective.dtlz.DTLZ3;
import org.uma.jmetal.problem.multiobjective.wfg.WFG1;
import org.uma.jmetal.runner.AbstractAlgorithmRunner;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.AlgorithmRunner;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.solutionattribute.impl.OverallConstraintViolation;

import com.mathworks.toolbox.javabuilder.MWException;

import MatlabPlot.LinePlot;
import MatlabPlot.ScatterPlot;

/**
 * Class for configuring and running the RSEA algorithm;
 * A Region Search Evolutionary Algorithm for Many-Objective Optimization
 *
 * @author Yongqi Liu <yongqil@hust.edu.cn>
 */
public class RSEARunner extends AbstractAlgorithmRunner {
	/**
	 * @param args Command line arguments.
	 * @throws SecurityException
	 * Invoking command:
  java org.uma.jmetal.runner.multiobjective.MOEADRunner problemName [referenceFront]
	 */

	private static int nfun = 3;

	public static void main(String[] args) throws FileNotFoundException, MWException {
		
		DoubleProblem problem;
		Algorithm<List<DoubleSolution>> algorithm;
		MutationOperator<DoubleSolution> mutation;
		CrossoverOperator<DoubleSolution> crossover;
		String referenceParetoFront = "";
		//      referenceParetoFront = "/pareto_fronts/DTLZ1.10D.pf";
		//      referenceParetoFront = "/pareto_fronts/Convex_DTLZ2.3D.pf";


//		problem = new  DTLZ1(nfun+4,nfun);
//		problem = new  DTLZ2(nfun+9,nfun);
//		problem = new  DTLZ3(nfun+9,nfun);
		//    problem = new  DTLZ4(nfun+9,nfun);
		//    problem = new  DTLZ5(nfun+9,nfun);
		//    problem = new  DTLZ6(nfun+9,nfun);
		//    problem = new  DTLZ7(nfun+19,nfun);
//		problem = new  Convex_DTLZ2(nfun+9,nfun);
		//    problem = new  Convex_DTLZ3(nfun+9,nfun);
		//    problem = new  SDTLZ1(nfun+4,nfun);
		//    problem = new  SDTLZ2(nfun+9,nfun);
		//    problem = new  ConvexC2_DTLZ2(nfun+9,nfun);
//		    problem = new  C1_DTLZ1(nfun+4,nfun);
//		    problem = new  C1_DTLZ3(nfun+9,nfun);
//		    problem = new  C2_DTLZ2(nfun+9,nfun);
//		    problem = new  ConvexC2_DTLZ2(nfun+9,nfun);
//		    problem = new  C3_DTLZ1(nfun+4,nfun,nfun);
//		    problem = new  C3_DTLZ4(nfun+4,nfun,nfun);
		    problem = new  WFG1(2*(nfun-1),20,nfun);
		//    problem = new  WFG2(2*(nfun-1),20,nfun);
		//    problem = new  WFG3(2*(nfun-1),20,nfun);
		//    problem = new  WFG4(2*(nfun-1),20,nfun);
		//    problem = new  WFG5(2*(nfun-1),20,nfun);
		//    problem = new  WFG6(2*(nfun-1),20,nfun);
		//    problem = new  WFG7(2*(nfun-1),20,nfun);
		//    problem = new  WFG8(2*(nfun-1),20,nfun);
		//    problem = new  WFG9(2*(nfun-1),20,nfun);

		double cr = 0.5 ;
		double f = 0.5 ;
		crossover = new DifferentialEvolutionCrossover(cr, f, "rand/1/bin");

		double mutationProbability = 1.0 / problem.getNumberOfVariables();
		double mutationDistributionIndex = 20.0;
		mutation = new PolynomialMutation(mutationProbability, mutationDistributionIndex);

		algorithm =new MOEADBuilder(problem, MOEADBuilder.Variant.RSEA)
				.setCrossover(crossover)
				.setMutation(mutation)
				.setPopulationSize(91)
				.setMaxEvaluations(91*3000)
				.setNeighborhoodSelectionProbability(0.9)
				.setNeighborSize(20)
				.setPlot(true) // The plot function need matlab 
				.setNormalize(true)
				.build() ;

		AlgorithmRunner algorithmRunner = new AlgorithmRunner.Executor(algorithm)
				.execute() ;
		List<DoubleSolution> population = algorithm.getResult() ;
		long computingTime = algorithmRunner.getComputingTime() ;
		JMetalLogger.logger.info("Total execution time: " + computingTime + "ms");

		printFinalSolutionSet(population);
		if (!referenceParetoFront.equals("")) {
			printQualityIndicators(population, referenceParetoFront) ;
		}
		if(false){
			OverallConstraintViolation<DoubleSolution> overallConstraintViolation = new OverallConstraintViolation<DoubleSolution>(); 
			ScatterPlot scatterplot = new ScatterPlot();
			LinePlot linePlotVar = new LinePlot();
			LinePlot linePlotObj = new LinePlot();
			double[][] fit = new double[population.get(0).getNumberOfObjectives()][population.size()];
			double[][] var = new double[population.size()][population.get(0).getNumberOfVariables()];
			for (int i = 0; i < population.size(); i++) {
				if(overallConstraintViolation.getAttribute(population.get(i))<0) continue;
				for (int j = 0; j < population.get(0).getNumberOfObjectives(); j++) {
					fit[j][i] = (population.get(i).getObjective(j));
				}
				for (int j = 0; j < population.get(0).getNumberOfVariables(); j++) {
					var[i][j] = (double) population.get(i).getVariableValue(j);
				}
			}
			if(fit.length==2){
				scatterplot.plot2scatter(fit[0],fit[1]);
			}else if(fit.length==3){

				scatterplot.plot3scatter(fit[0],fit[1],fit[2]);

			}else{
				linePlotObj.plotlines(fit);
			}
			linePlotVar.plotlines(var);
		}



	}
}
