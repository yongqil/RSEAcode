package org.mine.uitl;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.mine.exec.RunConfiguration;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.JMetalException;

public final class FileUtil
{
	/** 
	 * 读取文件并按行输出
	 * @param filePath
	 * @param spec 允许解析的最大行数， spec==null时，解析所有行
	 * @return
	 * @author l00428364
	 * @since 2017-12-8
	 */
	public static String[] read(final String filePath, final Integer spec)
	{
		File file = new File(filePath);
		// 当文件不存在或者不可读时
		if ((!isFileExists(file)) || (!file.canRead()))
		{
			System.out.println("file [" + filePath + "] is not exist or cannot read!!!");
			return null;
		}

		List<String> lines = new LinkedList<String>();
		BufferedReader br = null;
		FileReader fb = null;
		try
		{
			fb = new FileReader(file);
			br = new BufferedReader(fb);

			String str = null;
			int index = 0;
			while (((spec == null) || index++ < spec) && (str = br.readLine()) != null)
			{
				lines.add(str);
			}
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		finally
		{
			closeQuietly(br);
			closeQuietly(fb);
		}

		return lines.toArray(new String[lines.size()]);
	}
	/** 
	 * 写文件
	 * @param filePath 输出文件路径
	 * @param content 要写入的内容
	 * @param append 是否追加
	 * @return
	 * @author l00428364
	 * @since 2017-12-8
	 */
	public static int write(final String filePath, final String[] contents, final boolean append)
	{
		File file = new File(filePath);
		if (contents == null)
		{
			System.out.println("file [" + filePath + "] invalid!!!");
			return 0;
		}

		// 当文件存在但不可写时
		if (isFileExists(file) && (!file.canRead()))
		{
			return 0;
		}

		FileWriter fw = null;
		BufferedWriter bw = null;
		try
		{
			if (!isFileExists(file))
			{
				file.createNewFile();
			}

			fw = new FileWriter(file, append);
			bw = new BufferedWriter(fw);
			for (String content : contents)
			{
				if (content == null)
				{
					continue;
				}
				bw.write(content);
				bw.newLine();
			}
		}
		catch (IOException e)
		{
			e.printStackTrace();
			return 0;
		}
		finally
		{
			closeQuietly(bw);
			closeQuietly(fw);
		}

		return 1;
	}

	private static void closeQuietly(Closeable closeable)
	{
		try
		{
			if (closeable != null)
			{
				closeable.close();
			}
		}
		catch (IOException e)
		{
		}
	}

	private static boolean isFileExists(final File file)
	{
		if (file.exists() && file.isFile())
		{
			return true;
		}

		return false;
	}

	public List<RunConfiguration> readConfigurations(String filePath) {
		List<RunConfiguration> configurations = new ArrayList<>();
		InputStream inputStream = getClass().getResourceAsStream(filePath);
		try {
			InputStreamReader isr = new InputStreamReader(inputStream);
			BufferedReader br = new BufferedReader(isr);
			String lineTxt = null;
			while ((lineTxt = br.readLine()) != null) {
				String line[] = lineTxt.split("\t");
				RunConfiguration configuration = new RunConfiguration();
				configuration.setProblemName(line[0]);
				configuration.setFunNum(Integer.parseInt(line[1]));
				configuration.setPopsize(Integer.parseInt(line[2]));
				configuration.setRunNum(Integer.parseInt(line[3]));
				configuration.setNeighborhoodSelectionProbability(Double.parseDouble(line[4]));
				configuration.setNeighborhoodSize(Integer.parseInt(line[5]));
				configuration.setCalIGD(Boolean.parseBoolean(line[6]));
				configuration.setReferenceFrontsPath(line[7]);
				configuration.setNormalized(Boolean.parseBoolean(line[8]));
				configuration.newProblem();
				configurations.add(configuration);
				System.out.println(lineTxt);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return configurations;
	}

	public void writeFunction(String filePath, List<List<DoubleSolution>> populations) {
		try {
			File filename = new File(filePath);  
			if (!filename.exists()) {  
				filename.createNewFile();  
			} 
			FileOutputStream outputStream = new FileOutputStream(filePath);
			OutputStreamWriter outputStreamWriter = new OutputStreamWriter(outputStream);
			BufferedWriter bufferedWriter = new BufferedWriter(outputStreamWriter);

			for (int n = 0; n < populations.size(); n++) {
				bufferedWriter.write("#");
				bufferedWriter.newLine();
				List<DoubleSolution> solutionList = populations.get(n);
				if (solutionList.size() > 0) {
					int numberOfObjectives = solutionList.get(0).getNumberOfObjectives();
					for (int i = 0; i < solutionList.size(); i++) {
						for (int j = 0; j < numberOfObjectives; j++) {
							bufferedWriter.write(solutionList.get(i).getObjective(j) + "\t");
						}
						bufferedWriter.newLine();
					}
				}
			}
			bufferedWriter.write("#");
			bufferedWriter.close();
		} catch (IOException e) {
			throw new JMetalException("Error printing objecives to file: ", e);
		}

	}

	public void writeIGD(String filePath, List<Double> iGDs) {
		try {
			File filename = new File(filePath);  
			if (!filename.exists()) {  
				filename.createNewFile();  
			} 
			FileOutputStream outputStream = new FileOutputStream(filePath);
			OutputStreamWriter outputStreamWriter = new OutputStreamWriter(outputStream);
			BufferedWriter bufferedWriter = new BufferedWriter(outputStreamWriter);

			if (iGDs.size() > 0) {
				for (int i = 0; i < iGDs.size(); i++) {
					bufferedWriter.write(iGDs.get(i) + "");
					bufferedWriter.newLine();
				}
			}
			bufferedWriter.close();
		} catch (IOException e) {
			throw new JMetalException("Error printing objecives to file: ", e);
		}
	}

}
