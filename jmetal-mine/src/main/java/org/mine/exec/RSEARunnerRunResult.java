package org.mine.exec;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

import org.mine.uitl.FileUtil;
import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.algorithm.multiobjective.moead.AbstractMOEAD;
import org.uma.jmetal.algorithm.multiobjective.moead.MOEADBuilder;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.operator.impl.crossover.DifferentialEvolutionCrossover;
import org.uma.jmetal.operator.impl.mutation.PolynomialMutation;
import org.uma.jmetal.problem.DoubleProblem;
import org.uma.jmetal.qualityindicator.impl.InvertedGenerationalDistancePlus;
import org.uma.jmetal.runner.AbstractAlgorithmRunner;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.AlgorithmRunner;
import org.uma.jmetal.util.front.Front;
import org.uma.jmetal.util.front.imp.ArrayFront;
import org.uma.jmetal.util.front.util.FrontNormalizer;
import org.uma.jmetal.util.front.util.FrontUtils;
import org.uma.jmetal.util.point.util.PointSolution;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;
import org.uma.jmetal.util.pseudorandom.impl.JavaRandomGenerator;

public class RSEARunnerRunResult extends AbstractAlgorithmRunner {
	/**
	 * @param args Command line arguments.
	 * @throws SecurityException
	 * Invoking command:
  java org.uma.jmetal.runner.multiobjective.MOEADRunner problemName [referenceFront]
	 */


	public static void main(String[] args) throws FileNotFoundException {

		Algorithm<List<DoubleSolution>> algorithm;
		DoubleProblem problem;
		MutationOperator<DoubleSolution> mutation;
		CrossoverOperator<DoubleSolution> crossover;
		SelectionOperator<List<DoubleSolution>, DoubleSolution> selection;
		// read configration
		List<RunConfiguration> configurations = new FileUtil().readConfigurations("/runConfiguration.txt");
		JMetalRandom randomGenerator = JMetalRandom.getInstance() ;
		
		
//		randomGenerator.setSeed(1521209757553l);//

		for (RunConfiguration runConfiguration : configurations) {
			String referenceParetoFront = runConfiguration.getReferenceFrontsPath();
			int popsize = runConfiguration.getPopsize();
			int runNum = runConfiguration.getRunNum();
			problem  = runConfiguration.getProblem();
			
			
			double cr = 0.5 ;
		    double f = 0.5 ;
		    crossover = new DifferentialEvolutionCrossover(cr, f, "rand/1/bin");
		    

			double mutationProbability = 1.0 / problem.getNumberOfVariables();
			double mutationDistributionIndex = 20.0;
			mutation = new PolynomialMutation(mutationProbability, mutationDistributionIndex);

		    
			algorithm =new MOEADBuilder(problem, MOEADBuilder.Variant.RSEA)
					.setCrossover(crossover)
					.setMutation(mutation)
					.setPopulationSize(popsize)
					.setMaxEvaluations(popsize*runNum)
					.setResultPopulationSize(popsize)
					.setNeighborhoodSelectionProbability(runConfiguration.getNeighborhoodSelectionProbability())
					.setMaximumNumberOfReplacedSolutions(1)
					.setNeighborSize(runConfiguration.getNeighborhoodSize())
					.setFunctionType(AbstractMOEAD.FunctionType.PBI)
					.setDataDirectory("MOEAD_Weights")
					.setNormalize(runConfiguration.isNormalized())
					.build() ;
			
			List<List<DoubleSolution>> populations = new ArrayList<>();
			List<Double> IGDs = new ArrayList<>();
			for (int i = 0; i < 20; i++) {
				randomGenerator.setRandomGenerator( new JavaRandomGenerator() );
				System.out.println(randomGenerator.getSeed());
				AlgorithmRunner algorithmRunner = new AlgorithmRunner.Executor(algorithm)
						.execute() ;
				List<DoubleSolution> population = algorithm.getResult() ;
				populations.add(population);
				
//				long computingTime = algorithmRunner.getComputingTime() ;
//				JMetalLogger.logger.info("Total execution time: " + computingTime + "ms");
//				printFinalSolutionSet(population);
				System.out.println("Problem name :"+runConfiguration.getProblemName()+"\tNumber of objective:"+runConfiguration.getFunNum());
				if (runConfiguration.isCalIGD() && !referenceParetoFront.equals("")) {
					Front referenceFront = new ArrayFront(referenceParetoFront);
					FrontNormalizer frontNormalizer = new FrontNormalizer(referenceFront) ;

				    Front normalizedReferenceFront = frontNormalizer.normalize(referenceFront) ;
				    Front normalizedFront = frontNormalizer.normalize(new ArrayFront(population)) ;
				    List<PointSolution> normalizedPopulation = FrontUtils
				        .convertFrontToSolutionList(normalizedFront) ;
				    Double igd;
				    if(runConfiguration.getProblemName().startsWith("S")){
				    	igd = new InvertedGenerationalDistancePlus<PointSolution>(normalizedReferenceFront).evaluate(normalizedPopulation);
				    }else{
				    	igd = new InvertedGenerationalDistancePlus<DoubleSolution>(referenceFront).evaluate(population);
				    }
					IGDs.add(igd);
					System.out.println(igd);
				}
			}
			new FileUtil().writeFunction("result-RSEA/"+runConfiguration.getProblemName()+"_"+runConfiguration.getFunNum()+"_"
			+runConfiguration.getNeighborhoodSelectionProbability()+"_"+runConfiguration.getNeighborhoodSize()+"_FUN.txt",populations);
			if (runConfiguration.isCalIGD() && !referenceParetoFront.equals("")) {
				new FileUtil().writeIGD("result-RSEA/"+runConfiguration.getProblemName()+"_"+runConfiguration.getFunNum()+"_"
						+runConfiguration.getNeighborhoodSelectionProbability()+"_"+runConfiguration.getNeighborhoodSize()+"_IGD.txt",IGDs);
			}
			
		}

	}

	


}
