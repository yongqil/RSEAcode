package org.mine.uitl;
/**
 * @author 刘永琦 @qq 510596430 @tel 136168681951
 * @Email 510596430@qq.com
 * @Date 2016年8月1日 下午4:54:19
 * @version 1.0
 * @From 华中科技大学
 */
public interface SimpleFunctionMethod {

	/**
	 * 差值
	 * @param X 原序列X
	 * @param Y 原序列Y
	 * @param y 现需要求x,已知y
	 * @return 返回x
	 */
	double chazhi(double[] X, double[] Y, double y);

	double[][] GenerateReferencePoints(int M, int P);
	
}
