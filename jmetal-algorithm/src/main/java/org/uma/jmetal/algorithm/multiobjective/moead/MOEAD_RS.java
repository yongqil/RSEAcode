package org.uma.jmetal.algorithm.multiobjective.moead;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import org.uma.jmetal.algorithm.multiobjective.moead.util.MOEADUtils;
import org.uma.jmetal.algorithm.multiobjective.moead.util.SimpleFunctionMethod;
import org.uma.jmetal.algorithm.multiobjective.moead.util.SimpleFunctionMethodImpl;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.impl.crossover.DifferentialEvolutionCrossover;
import org.uma.jmetal.operator.impl.crossover.SBXCrossover;
import org.uma.jmetal.problem.ConstrainedProblem;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.comparator.DominanceComparator;
import org.uma.jmetal.util.comparator.impl.OverallConstraintViolationComparator;


/**
 * @author 刘永琦 @qq 510596430 @tel 136168681951
 * @Email 510596430@qq.com
 * @Date 2017年12月20日 上午8:44:41
 * @version 1.0
 * @param <S>
 * @From 华中科技大学
 */
@SuppressWarnings("serial")
public class MOEAD_RS<S extends Solution<?>>  extends AbstractMOEAD<DoubleSolution>{

	protected  CrossoverOperator<DoubleSolution>  crossover ;

	private static SimpleFunctionMethod Smethod = new SimpleFunctionMethodImpl();

	private static final Comparator<Solution<?>> DOMINANCE_COMPARATOR = new DominanceComparator<Solution<?>>();
	private static final Comparator<Solution<?>> CONSTRAINT_VIOLATION_COMPARATOR =
			new OverallConstraintViolationComparator<Solution<?>>();

	boolean normalnized;
	List<DoubleSolution> extremePoint = new ArrayList<>();
	List<Double> intercepts;
 	
	public MOEAD_RS(Problem<DoubleSolution> problem,
			int populationSize,
			int resultPopulationSize,
			int maxEvaluations,
			MutationOperator<DoubleSolution> mutation,
			CrossoverOperator<DoubleSolution> crossover,
			FunctionType functionType,
			String dataDirectory,
			double neighborhoodSelectionProbability,
			int maximumNumberOfReplacedSolutions,
			int neighborSize, boolean plot, boolean normalnized) {
		super(problem, populationSize, resultPopulationSize, maxEvaluations, crossover, mutation, functionType,
				dataDirectory, neighborhoodSelectionProbability, maximumNumberOfReplacedSolutions,
				neighborSize);
		this.plot=plot;
		this.crossover = crossover;
		this.normalnized=normalnized;
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

			int[] permutation = new int[populationSize];
			MOEADUtils.randomPermutation(permutation, populationSize);

			double F = 0,CR = 0;
			for (int i = 0; i < populationSize; i++) {
				int subProblemId = permutation[i];
				NeighborType neighborType = chooseNeighborType() ;
				List<DoubleSolution> children;
				List<DoubleSolution> parents =  new ArrayList<>(2);
				if(crossover instanceof DifferentialEvolutionCrossover){
					
					F = randomGenerator.nextDouble()<0.1? nextGaussian((double) population.get(i).getAttribute("F"),0.1) : (double) population.get(i).getAttribute("F");
					CR = randomGenerator.nextDouble()<0.1? nextGaussian((double) population.get(i).getAttribute("CR"),0.1) : (double) population.get(i).getAttribute("CR");
					if(F<=0 || F>1){
						F = (double) population.get(i).getAttribute("F");
					}
					if(CR<=0 || CR>1){
						CR = (double) population.get(i).getAttribute("CR");
					}
					((DifferentialEvolutionCrossover) crossover).setF(F);
					((DifferentialEvolutionCrossover) crossover).setCr(CR);
					parents = parentSelection(subProblemId, neighborType) ;
					((DifferentialEvolutionCrossover) crossover).setCurrentSolution(population.get(subProblemId));
				}else if(crossover instanceof SBXCrossover){
					List<Integer> matingPool= matingSelection(subProblemId, 2, neighborType);
					parents.add(population.get(matingPool.get(0)));
					parents.add(population.get(matingPool.get(1)));
				}

				children = crossover.execute(parents);
				DoubleSolution child = children.get(0) ;
				mutationOperator.execute(child);
				problem.evaluate(child);
				if (problem instanceof ConstrainedProblem) {
					((ConstrainedProblem<DoubleSolution>) problem).evaluateConstraints(child);
				}
				child.setAttribute("F", F);
				child.setAttribute("CR", CR);
				evaluations++;
				updateIdealPoint(child);
				if(normalnized){
					updateNadirpoint(child);
				}
				
				int childReg = findRegion(child);
				
				updateNeighborhoodNew(child, subProblemId, neighborType, childReg, 1);
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

	private double nextGaussian(double mean,double std_dev){
	    return mean+(nextGaussian()*std_dev);
	}

	private double nextGaussian(){
	    double u=0.0, v=0.0, w=0.0, c=0.0;
	    do{
	        u=randomGenerator.nextDouble()*2-1.0;
	        v=randomGenerator.nextDouble()*2-1.0;
	        w=u*u+v*v;
	    }while(w==0.0||w>=1.0);
	    c=Math.sqrt((-2*Math.log(w))/w);
	    //return [u*c,v*c];
	    return u*c;
	}
	
	private void updateNeighborhoodNew(DoubleSolution individual, int subProblemId,
			org.uma.jmetal.algorithm.multiobjective.moead.AbstractMOEAD.NeighborType neighborType, int childReg,
			int type) {
		int size;
		if (neighborType == NeighborType.NEIGHBOR) {
			size = neighborhood[subProblemId].length;
		} else {
			size = population.size();
		}
		int[] perm = new int[size];
		MOEADUtils.randomPermutation(perm, size);

		for (int i = 0; i < size; i++) {
			int k;
			if (neighborType == NeighborType.NEIGHBOR) {
				k = neighborhood[subProblemId][perm[i]];
			} else {
				k = perm[i];
			}
			int regionC;
			double cosP, cosC;
			double f1, f2;
			
			regionC = childReg;
			cosP = cosFunction2(population.get(k), lambda[k]);
			cosC = cosFunction2(individual, lambda[k]);
			double	flagDominate = CONSTRAINT_VIOLATION_COMPARATOR.compare( population.get(k), individual);
			if ( flagDominate==1 )   {
				population.set(k, (DoubleSolution)individual.copy());
			}
			
			if(regionC!=k && cosP>cosC){
					continue;
			}else if(type==1){
					f1 = fitFunction(population.get(k), lambda[k]);
					f2 = fitFunction(individual, lambda[k]);
					if(f2 < f1){
						population.set(k, (DoubleSolution)individual.copy());
					}
			}else if(type==2){
				if (flagDominate == 0) {
					flagDominate = DOMINANCE_COMPARATOR.compare( population.get(k), individual);
				}
				if ( flagDominate==1 )   {
					population.set(k, (DoubleSolution)individual.copy());
				}else if(flagDominate==0){
					f1 = calculateDistance2(population.get(k), lambda[k], idealPoint, nadirPoint);
					f2 = calculateDistance2(individual, lambda[k], idealPoint, nadirPoint);
					if(f2 < f1 ){
						population.set(k, (DoubleSolution)individual.copy());
					}
				}
			}
		}
	}



	// ----------------------------------------------------------------------
	// ASF: Achivement Scalarization Function
	// I implement here a effcient version of it, which only receives the index
	// of the objective which uses 1.0; the rest will use 0.00001. This is 
	// different to the one impelemented in C++
	// ----------------------------------------------------------------------
	private double ASF(DoubleSolution s, int index, double[] min) {
		double max_ratio = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < s.getNumberOfObjectives(); i++) {
			double weight = (index == i) ? 1.0 : 0.000001;
			max_ratio = Math.max(max_ratio, Math.abs((s.getObjective(i)-min[i]))/weight);
		}
		return max_ratio;
	}
	
	private List<DoubleSolution> findExtremePoints(List<DoubleSolution> population, double[] min) {
		List<DoubleSolution> extremePoints = new ArrayList<>();
		DoubleSolution min_indv = null;
		for (int f=0; f < problem.getNumberOfObjectives(); f+=1)
		{
			double min_ASF = Double.MAX_VALUE;	
			for (DoubleSolution s : population) { // only consider the individuals in the first front
				double asf = ASF(s, f, min);
				if ( asf < min_ASF ) {
					min_ASF = asf;
					min_indv = s;
				}
			}
			
			extremePoints.add(min_indv);
		}
		return extremePoints;
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
		boolean duplicate = false;

		List<Double> intercepts = new ArrayList<>();
		
		if (duplicate) // cannot construct the unique hyperplane (this is a casual method to deal with the condition)
		{
			List<Double> max_objs = FindMaxObjectives(population,child);
			for (int f=0; f<problem.getNumberOfObjectives(); f+=1)
			{
				// extreme_points[f] stands for the individual with the largest value of objective f
				intercepts.add(max_objs.get(f));
			}
		}else{
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
					return intercepts;
				}
				intercepts.add(1.0/x.get(f)+idealPoint[f]);
			}
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

	private int findRegion(DoubleSolution indiv) {
		double cos = -cosFunction2(indiv, lambda[0]);
		double maxCos=cos;
		int regionIdx = 0;
		for (int j = 1; j < populationSize; j++) {
			cos = -cosFunction2(indiv, lambda[j]);
			if(cos<maxCos){
				maxCos = cos;
				regionIdx = j;
			}
		}
		return regionIdx;
	}

	protected void initializePopulation() {
		population = new ArrayList<DoubleSolution>(populationSize);
		for (int i = 0; i < populationSize; i++) {
			DoubleSolution newSolution;
			newSolution = (DoubleSolution)problem.createSolution();
			newSolution.setAttribute("F", 0.5 /*randomGenerator.nextDouble(0, 1)*/);
			newSolution.setAttribute("CR",0.5);
			if (problem instanceof ConstrainedProblem) {
				 problem.evaluate(newSolution);
				((ConstrainedProblem<DoubleSolution>) problem).evaluateConstraints(newSolution);
				population.add(newSolution);
		      } else {
		        problem.evaluate(newSolution);
		        population.add(newSolution);
		      }
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
//		randomGenerator.setSeed(1527150823434l);//
		
	}


	/**
	 * theta = a/(sqrt(b)*sqrt(c))
	 */

	double cosFunction2(DoubleSolution individual, double[] lambda) throws JMetalException {
		
		double cos;
		double a=0,b=0,c=0,lmd,ob;
		for (int i = 0; i < lambda.length; i++) {
			lmd = lambda[i];
			if(normalnized){
				if(nadirPoint[i]-idealPoint[i] > 10e-10){
					ob = (individual.getObjective(i)-idealPoint[i])/(nadirPoint[i]-idealPoint[i]);
				}else{
					ob = (individual.getObjective(i)-idealPoint[i])/ 10e-10;
				}
			}else{
				ob = (individual.getObjective(i)-idealPoint[i]);
			}
			a += ob*lmd;
			b += ob*ob;
			c += lmd*lmd;
		}
		if(b==0){
			cos=1;
		}else{
			cos = a/(Math.sqrt(b)*Math.sqrt(c));
		}
		return cos;
	}
	

	  public double calculateDistance2(DoubleSolution indiv, double[] lambda,
	                                   double[] z_, double[] nz_) {
			double d1, d2, nl;
			d1 = d2 = nl = 0.0;

			for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
				if(normalnized){
					if(nadirPoint[i]-idealPoint[i] > 10e-10){
						d1 += (indiv.getObjective(i) - idealPoint[i])/(nadirPoint[i]-idealPoint[i]) * lambda[i];
					}else{
						d1 += (indiv.getObjective(i) - idealPoint[i])/10e-10 * lambda[i];
					}
				}else{
					d1 += (indiv.getObjective(i) - idealPoint[i]) * lambda[i];
				}
				nl += Math.pow(lambda[i], 2.0);
			}
			nl = Math.sqrt(nl);
			d1 = Math.abs(d1) / nl;

			for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
				if(normalnized){
					if(nadirPoint[i]-idealPoint[i] > 10e-10){
						d2 += Math.pow((indiv.getObjective(i) - idealPoint[i])/(nadirPoint[i]-idealPoint[i]) - d1 * (lambda[i] / nl), 2.0);
					}else{
						d2 += Math.pow((indiv.getObjective(i) - idealPoint[i])/10e-10 - d1 * (lambda[i] / nl), 2.0);
					}
				}else{
					d2 += Math.pow((indiv.getObjective(i) - idealPoint[i]) - d1 * (lambda[i] / nl), 2.0);
				}
			}
			d2 = Math.sqrt(d2);

			return d2;
	  }

	
	double fitFunction(DoubleSolution individual, double[] lambda) throws JMetalException {
		double fitness;
		double d1, d2, nl;
		double theta = 5.0;
		double diff = 0 ;
		d1 = d2 = nl = 0.0;

		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
			if(normalnized){
				if(nadirPoint[i]-idealPoint[i] > 10e-10){
					diff = (individual.getObjective(i) - idealPoint[i])/(nadirPoint[i]-idealPoint[i]) ;
				}else{
					diff = (individual.getObjective(i) - idealPoint[i])/ 10e-10 ;
				}
			}else{
				diff = (individual.getObjective(i) - idealPoint[i]);
			}
			d1 += diff * lambda[i];
			nl += Math.pow(lambda[i], 2.0);
		}
		nl = Math.sqrt(nl);
		d1 = Math.abs(d1) / nl;

		for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
			if(normalnized){
				if(nadirPoint[i]-idealPoint[i] > 10e-10){
					diff = (individual.getObjective(i) - idealPoint[i])/(nadirPoint[i]-idealPoint[i]);
				}else{
					diff = (individual.getObjective(i) - idealPoint[i])/ 10e-10;
				}
			}else{
				diff = (individual.getObjective(i) - idealPoint[i]);
			}
			d2 += Math.pow(diff - d1 * (lambda[i] / nl), 2.0);
		}
		d2 = Math.sqrt(d2);

		fitness = (d1 + theta * d2);
		return fitness;
	}

	@Override public String getName() {
		return "MOEAD_RS" ;
	}

	@Override public String getDescription() {
		return "Multi-Objective Evolutionary Algorithm based on Decomposition with region search strategy" ;
	}

	private boolean isPlot() {
		return plot;
	}
	
	
}
