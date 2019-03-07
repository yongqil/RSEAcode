package org.mine.exec;

import java.io.Serializable;

import org.uma.jmetal.problem.DoubleProblem;
import org.uma.jmetal.problem.multiobjective.dtlz.Convex_DTLZ2;
import org.uma.jmetal.problem.multiobjective.dtlz.DTLZ1;
import org.uma.jmetal.problem.multiobjective.dtlz.DTLZ2;
import org.uma.jmetal.problem.multiobjective.dtlz.DTLZ3;
import org.uma.jmetal.problem.multiobjective.dtlz.DTLZ4;
import org.uma.jmetal.problem.multiobjective.dtlz.SDTLZ1;
import org.uma.jmetal.problem.multiobjective.dtlz.SDTLZ2;
import org.uma.jmetal.problem.multiobjective.wfg.WFG1;
import org.uma.jmetal.problem.multiobjective.wfg.WFG2;
import org.uma.jmetal.problem.multiobjective.wfg.WFG3;
import org.uma.jmetal.problem.multiobjective.wfg.WFG4;
import org.uma.jmetal.problem.multiobjective.wfg.WFG5;
import org.uma.jmetal.problem.multiobjective.wfg.WFG6;
import org.uma.jmetal.problem.multiobjective.wfg.WFG7;
import org.uma.jmetal.problem.multiobjective.wfg.WFG8;
import org.uma.jmetal.problem.multiobjective.wfg.WFG9;

/**
 * @author 刘永琦 @qq 510596430 @tel 136168681951
 * @Email 510596430@qq.com
 * @Date 2018年3月16日 下午4:10:42
 * @version 1.0
 * @From 华中科技大学
 */
public class RunConfiguration implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -749461499715662715L;

	private String problemName;

	private DoubleProblem problem;

	private int funNum;

	private int popsize;

	private int runNum;

	private double neighborhoodSelectionProbability;

	private int neighborhoodSize;

	private String referenceFrontsPath;
	
	private boolean normalized;

	private boolean calIGD;

	public DoubleProblem getProblem() {
		return problem;
	}

	public int getFunNum() {
		return funNum;
	}

	public int getPopsize() {
		return popsize;
	}

	public int getRunNum() {
		return runNum;
	}

	public double getNeighborhoodSelectionProbability() {
		return neighborhoodSelectionProbability;
	}

	public int getNeighborhoodSize() {
		return neighborhoodSize;
	}

	public void setProblem(DoubleProblem problem) {
		this.problem = problem;
	}

	public void setFunNum(int funNum) {
		this.funNum = funNum;
	}

	public void setPopsize(int popsize) {
		this.popsize = popsize;
	}

	public void setRunNum(int runNum) {
		this.runNum = runNum;
	}

	public void setNeighborhoodSelectionProbability(double neighborhoodSelectionProbability) {
		this.neighborhoodSelectionProbability = neighborhoodSelectionProbability;
	}

	public void setNeighborhoodSize(int neighborhoodSize) {
		this.neighborhoodSize = neighborhoodSize;
	}

	public boolean isCalIGD() {
		return calIGD;
	}

	public void setCalIGD(boolean calIGD) {
		this.calIGD = calIGD;
	}

	public String getReferenceFrontsPath() {
		return referenceFrontsPath;
	}

	public void setReferenceFrontsPath(String referenceFrontsPath) {
		this.referenceFrontsPath = referenceFrontsPath;
	}

	public String getProblemName() {
		return problemName;
	}

	public void setProblemName(String problemName) {
		this.problemName = problemName;
	}

	public void newProblem() {
		switch (problemName) {
		case "DTLZ1":
			problem = new DTLZ1(funNum+4,funNum);
			break;
		case "DTLZ2":
			problem = new DTLZ2(funNum+9,funNum);
			break;
		case "DTLZ3":
			problem = new DTLZ3(funNum+9,funNum);
			break;
		case "DTLZ4":
			problem = new DTLZ4(funNum+9,funNum);
			break;
		case "SDTLZ1":
			problem = new SDTLZ1(funNum+4,funNum);
			break;
		case "SDTLZ2":
			problem = new SDTLZ2(funNum+9,funNum);
			break;
		case "Convex_DTLZ2":
			problem = new Convex_DTLZ2(funNum+9,funNum);
			break;
		case "WFG1":
			problem = new WFG1(2*(funNum-1),20,funNum);;
			break;
		case "WFG2":
			problem = new WFG2(2*(funNum-1),20,funNum);;
			break;
		case "WFG3":
			problem = new WFG3(2*(funNum-1),20,funNum);;
			break;
		case "WFG4":
			problem = new WFG4(2*(funNum-1),20,funNum);;
			break;
		case "WFG5":
			problem = new WFG5(2*(funNum-1),20,funNum);;
			break;
		case "WFG6":
			problem = new WFG6(2*(funNum-1),20,funNum);;
			break;
		case "WFG7":
			problem = new WFG7(2*(funNum-1),20,funNum);;
			break;
		case "WFG8":
			problem = new WFG8(2*(funNum-1),20,funNum);;
			break;
		case "WFG9":
			problem = new WFG9(2*(funNum-1),20,funNum);;
			break;

		}

	}

	public boolean isNormalized() {
		return normalized;
	}

	public void setNormalized(boolean normalized) {
		this.normalized = normalized;
	}



}
