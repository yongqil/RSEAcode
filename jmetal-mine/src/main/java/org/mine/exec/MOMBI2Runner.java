package org.mine.exec;

import org.mine.uitl.FileUtil;
import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.algorithm.multiobjective.moead.AbstractMOEAD;
import org.uma.jmetal.algorithm.multiobjective.moead.MOEADBuilder;
import org.uma.jmetal.algorithm.multiobjective.mombi.MOMBI2;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.operator.impl.crossover.DifferentialEvolutionCrossover;
import org.uma.jmetal.operator.impl.crossover.SBXCrossover;
import org.uma.jmetal.operator.impl.mutation.PolynomialMutation;
import org.uma.jmetal.operator.impl.selection.BinaryTournamentSelection;
import org.uma.jmetal.problem.DoubleProblem;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.qualityindicator.impl.InvertedGenerationalDistancePlus;
import org.uma.jmetal.runner.AbstractAlgorithmRunner;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.AlgorithmRunner;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.ProblemUtils;
import org.uma.jmetal.util.comparator.RankingAndCrowdingDistanceComparator;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;
import org.uma.jmetal.util.front.Front;
import org.uma.jmetal.util.front.imp.ArrayFront;
import org.uma.jmetal.util.front.util.FrontNormalizer;
import org.uma.jmetal.util.front.util.FrontUtils;
import org.uma.jmetal.util.point.util.PointSolution;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;
import org.uma.jmetal.util.pseudorandom.impl.JavaRandomGenerator;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

/**
 * Class to configure and run the MOMBI2 algorithm
 *
 * @author Juan J. Durillo <juan@dps.uibk.ac.at>
 *
 * Reference: Improved Metaheuristic Based on the R2 Indicator for Many-Objective Optimization.
 * R. Hernández Gómez, C.A. Coello Coello. Proceeding GECCO '15 Proceedings of the 2015 on Genetic
 * and Evolutionary Computation Conference. Pages 679-686
 * DOI: 10.1145/2739480.2754776
 */
public class MOMBI2Runner extends AbstractAlgorithmRunner {
  /**
   * @param args Command line arguments.
   * @throws JMetalException
   * @throws FileNotFoundException
   * Invoking command:
    java org.uma.jmetal.runner.multiobjective.MOMBIRunner problemName [referenceFront]
   */
  public static void main(String[] args) throws JMetalException, FileNotFoundException {
	  
		Algorithm<List<DoubleSolution>> algorithm;
		DoubleProblem problem;
		MutationOperator<DoubleSolution> mutation;
		CrossoverOperator<DoubleSolution> crossover;
		SelectionOperator<List<DoubleSolution>, DoubleSolution> selection;
		// read configration
		List<RunConfiguration> configurations = new FileUtil().readConfigurations("/runConfiguration.txt");
		JMetalRandom randomGenerator = JMetalRandom.getInstance() ;
		
		String[] weights = {"weight_03D_12.sld","weight_05D.sld","weight_08D.sld","weight_10D.sld","weight_15D.sld"};
		
//		randomGenerator.setSeed(1521209757553l);//

		for (RunConfiguration runConfiguration : configurations) {
			String referenceParetoFront = runConfiguration.getReferenceFrontsPath();
			int popsize = runConfiguration.getPopsize();
			int runNum = runConfiguration.getRunNum();
			problem  = runConfiguration.getProblem();
			
			double crossoverProbability = 1.0 ;
			double crossoverDistributionIndex = 30.0 ;
			crossover = new SBXCrossover(crossoverProbability, crossoverDistributionIndex) ;
			
//			double cr = 0.5 ;
//		    double f = 0.5 ;
//		    crossover = new DifferentialEvolutionCrossover(cr, f, "rand/1/bin");
		    

			double mutationProbability = 1.0 / problem.getNumberOfVariables();
			double mutationDistributionIndex = 20.0;
			mutation = new PolynomialMutation(mutationProbability, mutationDistributionIndex);

	/*	    selection = new BinaryTournamentSelection<DoubleSolution>();
		    
		    algorithm = new NSGAIIIBuilder<DoubleSolution>(problem)
		    		.setPopulationSize(popsize)
		            .setCrossoverOperator(crossover)
		            .setMutationOperator(mutation)
		            .setSelectionOperator(selection)
		            .setMaxIterations(runNum)
		            .setPlot(false)
		            .build() ;
*/
		    
/*			algorithm =new MOEADBuilder(problem, MOEADBuilder.Variant.MOEADMyNew)
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
					.setPlot(false)
					.setNormalize(runConfiguration.isNormalized())
//					.setNormalize(true)
					.build() ;*/
			
			selection = new BinaryTournamentSelection<DoubleSolution>(new RankingAndCrowdingDistanceComparator<DoubleSolution>());

			int wIndex = 0;
			if(problem.getNumberOfObjectives()==3){
				wIndex=0;
			}else if(problem.getNumberOfObjectives()==5){
				wIndex=1;
			}else if(problem.getNumberOfObjectives()==8){
				wIndex=2;
			}else if(problem.getNumberOfObjectives()==10){
				wIndex=3;
			}else if(problem.getNumberOfObjectives()==15){
				wIndex=4;
			}
		    algorithm = new MOMBI2<>(problem,runNum,crossover,mutation,selection,new SequentialSolutionListEvaluator<DoubleSolution>(),
		    		"mombi2-weights/weight/"+weights[wIndex]);
		    
			List<List<DoubleSolution>> populations = new ArrayList<>();
			List<Double> IGDs = new ArrayList<>();
			for (int i = 0; i < 1; i++) {
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
				System.out.println(algorithmRunner.getComputingTime());
//				if (runConfiguration.isCalIGD() && !referenceParetoFront.equals("")) {
//					Front referenceFront = new ArrayFront(referenceParetoFront);
//					FrontNormalizer frontNormalizer = new FrontNormalizer(referenceFront) ;
//
//				    Front normalizedReferenceFront = frontNormalizer.normalize(referenceFront) ;
//				    Front normalizedFront = frontNormalizer.normalize(new ArrayFront(population)) ;
//				    List<PointSolution> normalizedPopulation = FrontUtils
//				        .convertFrontToSolutionList(normalizedFront) ;
//				    Double igd;
//				    if(runConfiguration.getProblemName().startsWith("S")){
//				    	igd = new InvertedGenerationalDistancePlus<PointSolution>(normalizedReferenceFront).evaluate(normalizedPopulation);
//				    }else{
//				    	igd = new InvertedGenerationalDistancePlus<DoubleSolution>(referenceFront).evaluate(population);
//				    }
//					IGDs.add(igd);
//					System.out.println(igd);
//				}
			}
//			new FileUtil().writeFunction("result-MOMIBI2/"+runConfiguration.getProblemName()+"_"+runConfiguration.getFunNum()+"_"
//			+runConfiguration.getNeighborhoodSelectionProbability()+"_"+runConfiguration.getNeighborhoodSize()+"_FUN.txt",populations);
//			if (runConfiguration.isCalIGD() && !referenceParetoFront.equals("")) {
//				new FileUtil().writeIGD("result-MOMIBI2/"+runConfiguration.getProblemName()+"_"+runConfiguration.getFunNum()+"_"
//						+runConfiguration.getNeighborhoodSelectionProbability()+"_"+runConfiguration.getNeighborhoodSize()+"_IGD.txt",IGDs);
//			}
			
		}
	  
  }
}
