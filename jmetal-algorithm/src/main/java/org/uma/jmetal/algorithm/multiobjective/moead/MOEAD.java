package org.uma.jmetal.algorithm.multiobjective.moead;

import java.util.ArrayList;
import java.util.List;

import org.uma.jmetal.algorithm.multiobjective.moead.util.MOEADUtils;
import org.uma.jmetal.algorithm.multiobjective.moead.util.SimpleFunctionMethod;
import org.uma.jmetal.algorithm.multiobjective.moead.util.SimpleFunctionMethodImpl;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.impl.crossover.DifferentialEvolutionCrossover;
import org.uma.jmetal.operator.impl.crossover.SBXCrossover;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.DoubleSolution;

import com.mathworks.toolbox.javabuilder.MWException;

import MatlabPlot.LinePlot;
import MatlabPlot.ScatterPlot;

/**
 * Class implementing the MOEA/D-DE algorithm described in :
 * Hui Li; Qingfu Zhang, "Multiobjective Optimization Problems With Complicated Pareto Sets,
 * MOEA/D and NSGA-II," Evolutionary Computation, IEEE Transactions on , vol.13, no.2, pp.284,302,
 * April 2009. doi: 10.1109/TEVC.2008.925798
 *
 * @author Antonio J. Nebro
 * @version 1.0
 */
@SuppressWarnings("serial")
public class MOEAD extends AbstractMOEAD<DoubleSolution> {
	protected  CrossoverOperator<DoubleSolution>  crossover ;
	protected ScatterPlot scatterplot = null;
	protected LinePlot linePlotVar = null;
	protected LinePlot linePlotObj = null;
	private static SimpleFunctionMethod Smethod = new SimpleFunctionMethodImpl();
	
	boolean normalnized;
	List<DoubleSolution> extremePoint = new ArrayList<>();
	List<Double> intercepts;
	
  public MOEAD(Problem<DoubleSolution> problem,
      int populationSize,
      int resultPopulationSize,
      int maxEvaluations,
      MutationOperator<DoubleSolution> mutation,
      CrossoverOperator<DoubleSolution> crossover,
      FunctionType functionType,
      String dataDirectory,
      double neighborhoodSelectionProbability,
      int maximumNumberOfReplacedSolutions,
      int neighborSize,
      boolean normalnized) {
    super(problem, populationSize, resultPopulationSize, maxEvaluations, crossover, mutation, functionType,
        dataDirectory, neighborhoodSelectionProbability, maximumNumberOfReplacedSolutions,
        neighborSize);

    this.crossover = crossover;
    this.normalnized=normalnized;
    super.setNormalnized(normalnized);
  }

  @Override public void run() {
	  initializeUniformWeight();
    initializePopulation() ;
    initializeNeighborhood();
    initializeIdealPoint() ;
    if(normalnized){
    	extremePoint = 		findExtremePoints(population,idealPoint);
		intercepts =constructHyperplane(population, extremePoint,population.get(0));
		for (int i = 0; i < nadirPoint.length; i++) {
			nadirPoint[i] = intercepts.get(i);
		}
    }

    evaluations = populationSize ;
    do {
    	/*Plot*/
    	try {
    		if(false){
    			if(scatterplot==null){
    				scatterplot = new ScatterPlot();
    				linePlotVar = new LinePlot();
    				linePlotObj = new LinePlot();
    			}
    			double[][] fit = new double[population.get(0).getNumberOfObjectives()][population.size()];
    			double[][] var = new double[population.size()][population.get(0).getNumberOfVariables()];
    			for (int i = 0; i < population.size(); i++) {
    				for (int j = 0; j < population.get(0).getNumberOfObjectives(); j++) {
    					fit[j][i] = (population.get(i).getObjective(j));
    				}
    				for (int j = 0; j < population.get(0).getNumberOfVariables(); j++) {
    					var[i][j] = (double) population.get(i).getVariableValue(j);
    				}
    			}
    			if(fit.length==2){
    				scatterplot.plot2scatter(fit[0],fit[1],"Gener"+evaluations/populationSize);
    			}else if(fit.length==3){

    				scatterplot.plot3scatter(fit[0],fit[1],fit[2],"Gener"+evaluations/populationSize);

    			}else{
    				linePlotObj.plotlines(fit,"Gener"+evaluations/populationSize);
    			}
    			linePlotVar.plotlines(var,"Gener"+evaluations/populationSize);
    		}
    	} catch (MWException e) {
    		// TODO Auto-generated catch block
    		e.printStackTrace();
    	}
    	
      int[] permutation = new int[populationSize];
      MOEADUtils.randomPermutation(permutation, populationSize);

      for (int i = 0; i < populationSize; i++) {
        int subProblemId = permutation[i];

        NeighborType neighborType = chooseNeighborType() ;
        List<DoubleSolution> parents =  new ArrayList<>(2);
		if(crossover instanceof DifferentialEvolutionCrossover){
			
			parents = parentSelection(subProblemId, neighborType) ;
			((DifferentialEvolutionCrossover) crossover).setCurrentSolution(population.get(subProblemId));
		}else if(crossover instanceof SBXCrossover){
			List<Integer> matingPool= matingSelection(subProblemId, 2, neighborType);
			parents.add(population.get(matingPool.get(0)));
			parents.add(population.get(matingPool.get(1)));
		}
        List<DoubleSolution> children = crossover.execute(parents);

        DoubleSolution child = children.get(0) ;
        mutationOperator.execute(child);
        problem.evaluate(child);

        evaluations++;

        updateIdealPoint(child);
        if(normalnized){
        	updateNadirpoint(child);
        }
        updateNeighborhood(child, subProblemId, neighborType);
      }
    } while (evaluations < maxEvaluations);

  }

  private void updateNadirpoint(DoubleSolution child) {
		boolean needUpdate=false;
		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
			if(ASF(child, i, idealPoint) < ASF(extremePoint.get(i), i, idealPoint)){
				extremePoint.set(i, (DoubleSolution) child.copy());
				needUpdate = true;
			}
		}
		if(needUpdate){
			intercepts =constructHyperplane(population, extremePoint,child);
			for (int i = 0; i < nadirPoint.length; i++) {
				nadirPoint[i] = intercepts.get(i);
			}
		}
		
	}
  private List<Double> FindMaxObjectives(List<DoubleSolution> population, DoubleSolution child) {
		List<Double> max_point = new ArrayList<>();
		for (int f=0; f<problem.getNumberOfObjectives(); f++){
			double max = child.getObjective(f);
			for (int i=0; i<population.size(); i++){
				if(population.get(i).getObjective(f) > max){
					max = population.get(i).getObjective(f);
				}
			}
			max_point.add(max);
		}

		return max_point;
	}
	
	public List<Double> constructHyperplane(List<DoubleSolution> population, List<DoubleSolution> extreme_points, DoubleSolution child) {
		List<Double> intercepts = new ArrayList<>();
		// Find the equation of the hyperplane
		List<Double> b = new ArrayList<>(); //(pop[0].objs().size(), 1.0);
		for (int i =0; i < problem.getNumberOfObjectives();i++)
			b.add(1.0);
		List<List<Double>> A=new ArrayList<>();
		for (DoubleSolution s : extreme_points)
		{
			List<Double> aux = new ArrayList<>();
			for (int i = 0; i < problem.getNumberOfObjectives(); i++)
				aux.add(s.getObjective(i)-idealPoint[i]);
			A.add(aux);
		}
		List<Double> x = guassianElimination(A, b);
		// Find intercepts
		for (int f=0; f<problem.getNumberOfObjectives(); f+=1)
		{
			if(x.get(f)<=0 || x.get(f).isNaN()){
				intercepts.clear();
				List<Double> max_objs = FindMaxObjectives(population,child);
				for (int f1=0; f1<problem.getNumberOfObjectives(); f1+=1)
				{
					intercepts.add(max_objs.get(f1));
				}
				/*for (int o=0; o<problem.getNumberOfObjectives(); o+=1)
				{
					intercepts.add(extreme_points.get(o).getObjective(o));
				}*/
				return intercepts;
			}
			intercepts.add(1.0/x.get(f)+idealPoint[f]);
		}
		return intercepts;
	}
	
	public List<Double> guassianElimination(List<List<Double>> A, List<Double> b) {
		List<Double> x = new ArrayList<>();

	    int N = A.size();
	    for (int i=0; i<N; i+=1)
	    {
	    	A.get(i).add(b.get(i));
	    }

	    for (int base=0; base<N-1; base+=1)
	    {
	        for (int target=base+1; target<N; target+=1)
	        {
	            double ratio = A.get(target).get(base)/A.get(base).get(base);
	            for (int term=0; term<A.get(base).size(); term+=1)
	            {
	                A.get(target).set(term, A.get(target).get(term) - A.get(base).get(term)*ratio);
	            }
	        }
	    }

	    for (int i = 0; i < N; i++)
	    	x.add(0.0);
	    
	    for (int i=N-1; i>=0; i-=1)
	    {
	        for (int known=i+1; known<N; known+=1)
	        {
	            A.get(i).set(N, A.get(i).get(N) - A.get(i).get(known)*x.get(known));
	        }
	        x.set(i, A.get(i).get(N)/A.get(i).get(i));
	    }
		return x;
	}

	private List<DoubleSolution> findExtremePoints(List<DoubleSolution> population, double[] min) {
		List<DoubleSolution> extremePoints = new ArrayList<>();
		DoubleSolution min_indv = null;
		for (int f=0; f < problem.getNumberOfObjectives(); f+=1)
		{
			double min_ASF = Double.MAX_VALUE;	
			for (DoubleSolution s : population) { // only consider the individuals in the first front
				double asf = ASF(s, f, min);
//				double asf = -COS(s, f, min);
//				double asf = D2(s, f, min);
				if ( asf < min_ASF ) {
					min_ASF = asf;
					min_indv = s;
				}
			}
			
			extremePoints.add(min_indv);
		}
		return extremePoints;
	}
	
	private double ASF(DoubleSolution s, int index, double[] min) {
		double max_ratio = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < s.getNumberOfObjectives(); i++) {
			double weight = (index == i) ? 1.0 : 0.000001;
			max_ratio = Math.max(max_ratio, Math.abs((s.getObjective(i)-min[i]))/weight);
		}
		return max_ratio;
	}

protected void initializePopulation() {
    population = new ArrayList<>(populationSize);
    for (int i = 0; i < populationSize; i++) {
      DoubleSolution newSolution = (DoubleSolution)problem.createSolution();

      problem.evaluate(newSolution);
      population.add(newSolution);
    }
  }
  
	@Override
	protected void initializeUniformWeight() {
		int  nfunc = problem.getNumberOfObjectives();
		this.lambda = new double[populationSize][nfunc];
		if(nfunc==2){
			for (int i = 0; i < populationSize; i++) {
				lambda[i][0] = (double)i/(populationSize-1) ==0  ? 0.0 : (double)i/(populationSize-1);
				lambda[i][1] = 1-(double)i/(populationSize-1) ==0  ? 0.0 : 1-(double)i/(populationSize-1);
			}
		}else {
			int p=12;
			if (nfunc==3){
				p=12;
				lambda=Smethod.GenerateReferencePoints(nfunc, p);
			}
			if (nfunc==4){
				p=7;
				lambda=Smethod.GenerateReferencePoints(nfunc, p);
			}
			if (nfunc==5){
				p=6;
				lambda=Smethod.GenerateReferencePoints(nfunc, p);
			}
			/*if (nfunc>=8 && nfunc<=10 ){
				int s=0;
				if(nfunc==8) s= 0;
				if(nfunc==10) s= 1;
				if(nfunc==15) s= 2;
				ExcelFactory factory = new ExcelFactory();
				*//**-----------------------读取------------------------------**//*
				Workbook wk = factory.readExcel("/referencepoint/w5~15.xls");
				Sheet  sheet = wk.getSheet(s);
				lambda = new double[sheet.getRows()][nfunc];
				for (int i = 0; i < sheet.getRows(); i++) {
					Cell[] cell = sheet.getRow(i); //把第一行返回成一个一个格子的形式给你
					for (int j = 0; j < cell.length; j++) {
						lambda[i][j] = Double.parseDouble(cell[j].getContents());
					}
				}
			}*/
			if(nfunc==8 || nfunc==10){
				double[][] lambda1 = Smethod.GenerateReferencePoints(nfunc, 3);
				double[][] lambda2 = Smethod.GenerateReferencePoints(nfunc, 2);
				for (int k = 0; k < lambda2.length; k++) {
					for (int j = 0; j < nfunc; j++) {
						lambda2[k][j] = 0.5/nfunc +0.5*lambda2[k][j];
					}
				}
				lambda = Smethod.concat(lambda1, lambda2);
				
			}
			if(nfunc==15){
				double[][] lambda1 = Smethod.GenerateReferencePoints(nfunc, 2);
				double[][] lambda2 = Smethod.GenerateReferencePoints(nfunc, 1);
				for (int k = 0; k < lambda2.length; k++) {
					for (int j = 0; j < lambda2.length; j++) {
						lambda2[k][j] = 0.5/nfunc +0.5*lambda2[k][j];
					}
				}
				lambda = Smethod.concat(lambda1, lambda2);
				
			}

		}
		populationSize = lambda.length;
		resultPopulationSize = populationSize;
		neighborhood = new int[populationSize][neighborSize];
//		randomGenerator.setSeed(1521097873054l);//
//		randomGenerator.setSeed(1521201544021l);//
//		randomGenerator.setSeed(1526369872320l);//
		
	}

  @Override public String getName() {
    return "MOEAD" ;
  }

  @Override public String getDescription() {
    return "Multi-Objective Evolutionary Algorithm based on Decomposition" ;
  }
}
