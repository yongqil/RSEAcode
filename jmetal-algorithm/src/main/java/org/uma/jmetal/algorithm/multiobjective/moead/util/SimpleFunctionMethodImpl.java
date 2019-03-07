package org.uma.jmetal.algorithm.multiobjective.moead.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author 刘永琦 @qq 510596430 @tel 136168681951
 * @Email 510596430@qq.com
 * @Date 2017年12月20日 上午9:02:22
 * @version 1.0
 * @From 华中科技大学
 */
public class SimpleFunctionMethodImpl implements SimpleFunctionMethod {

	/**
	 * 差值
	 * @param X 原序列X
	 * @param Y 原序列Y
	 * @param y 现需要求x,已知y
	 * @return 返回x
	 */
	public double chazhi(double[] X,double[] Y,double y) {
		int up=X.length-1,down=0,mid;
		double x;
		while(up-down>1){
			mid=Math.round((up+down)/2);
			if(y==Y[mid]) return x=X[mid];
			else if(y<=Y[mid]) up=mid;
			else down=mid;
		}
		x=X[down]+(y-Y[down])*(X[up]-X[down])/(Y[up]-Y[down]);	
		return x;
	}

	public double[][] GenerateReferencePoints(int M, int P) {
		List<List<Double>> ref = GetFixedRowSumIntegerMatrix(M,P);
		double[][] refPoint = new double[ref.size()][M];
		for (int i = 0; i < refPoint.length; i++) {
			for (int j = 0; j < refPoint[i].length; j++) {
				refPoint[i][j] = ref.get(i).get(j)/P;
			}
		}
		return refPoint;
	}

	private List<List<Double>> GetFixedRowSumIntegerMatrix(int M, int RowSum) {
		if (M == 1){
			Double a = (double) RowSum;
			List<List<Double>> refA = new ArrayList<List<Double>>();
			List<Double> A = new ArrayList<Double>();
			A.add(a);
			refA.add(A);
			return refA;
		}
		List<List<Double>> refA = new ArrayList<List<Double>>();
		for (int i = 0; i < RowSum+1; i++) {
			List<List<Double>> refB =GetFixedRowSumIntegerMatrix(M-1,RowSum-i);

			for (List<Double> b : refB) {
				b.add((double) i);
				refA.add(b);
			}
		}

		return refA;
	}

	@Override
	public <T> T[] concat(T[] first, T[] second) {
		T[] result = Arrays.copyOf(first, first.length + second.length);  
		System.arraycopy(second, 0, result, first.length, second.length);  
		return result;  
	}



}
