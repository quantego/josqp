package com.quantego.josqp;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.LinkedList;
import java.util.List;

public class Utils {
	
	public static double[] toDoubleArray(List<Double> array) {
        double[] result = new double[array.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = array.get(i);
        }
        return result;
    }

	public static double[] toDoubleArrayNeg(List<Double> array) {
		double[] result = new double[array.size()];
		for (int i = 0; i < result.length; i++) {
			result[i] = -array.get(i);
		}
		return result;
	}
	    
    public  static int[] toIntArray(List<Integer> array) {
        int[] result = new int[array.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = array.get(i);
        }
        return result;
    }
    
    public static double[] readDoubleColumn(String filename, String separator, int rowskip, int colindex) {
		List<Double> list = new LinkedList<>();
		try {
			FileInputStream in = new FileInputStream(filename);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			int k=0;
			while (k++<rowskip)
				br.readLine();
			String s;
			while ((s = br.readLine())!=null) {
				String[] sa = s.split(separator);
				if (!sa[colindex].isEmpty()) {
					if (sa[colindex].matches("-inf"))
						list.add(Double.NEGATIVE_INFINITY);
					else if (sa[colindex].matches("inf"))
						list.add(Double.POSITIVE_INFINITY);
					else
						list.add(Double.parseDouble(sa[colindex]));
				}
				else
					list.add(Double.NaN);
			}
			br.close();
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return toDoubleArray(list);
	}
    
    public static int[] readIntColumn(String filename, String separator, int rowskip, int colindex) {
		List<Integer> list = new LinkedList<>();
		try {
			FileInputStream in = new FileInputStream(filename);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			int k=0;
			while (k++<rowskip)
				br.readLine();
			String s;
			while ((s = br.readLine())!=null) {
				String[] sa = s.split(separator);
				if (!sa[colindex].isEmpty())
					list.add(Integer.parseInt(sa[colindex]));
				else
					list.add(Integer.MAX_VALUE);
			}
			br.close();
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return toIntArray(list);
	}

}
