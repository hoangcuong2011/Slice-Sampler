/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

import java.util.ArrayList;
import java.util.Random;

/**
 *
 * @author hoangcuong2011
 */
public class SliceSampler {
    class StatisticsPackage {

        ArrayList<Integer> cumusCounts = new ArrayList<>();
        ArrayList<Double> cumusProbs = new ArrayList<>();
        ArrayList<Double> cumusPoints = new ArrayList<>();
    }

    public StatisticsPackage cumuculativedistribution(ArrayList<Double> sample) {

        ArrayList<Integer> cumusCounts = new ArrayList<>();
        ArrayList<Double> cumusProbs = new ArrayList<>();
        ArrayList<Double> cumusPoints = new ArrayList<>();

        StatisticsPackage p = new StatisticsPackage();
        double min = 10000;
        for (int i = 0; i < sample.size(); i++) {
            if (sample.get(i) <= min) {
                min = sample.get(i);
            }
        }

        double max = -10000;
        for (int i = 0; i < sample.size(); i++) {
            if (sample.get(i) >= max) {
                max = sample.get(i);
            }
        }

        //pivot = max - min / 10;
        int bin_size = 20;
        double gaps = (max - min) / (double) bin_size;

        for (int b = 0; b < bin_size; b++) {
            int count = 0;
            double threshold = (min + gaps * (b + 1));
            cumusPoints.add(threshold);
            for (int i = 0; i < sample.size(); i++) {
                //if (sample.get(i) <= (threshold)) - we can lost precision here.
                if (threshold - sample.get(i) > -0.0001) {
                    count++;
                }
            }
            cumusCounts.add(count);
            cumusProbs.add((double) count / (double) sample.size());
        }
        /*for (int b = 0; b < bin_size; b++) {
            System.out.println(cumusCounts.get(b));
        }*/
        p.cumusCounts = cumusCounts;
        p.cumusPoints = cumusPoints;
        p.cumusProbs = cumusProbs;
        return p;
    }

    public StatisticsPackage cumuculativedistribution(ArrayList<Double> sample, ArrayList<Double> cumusPoints_for_reference) {

        ArrayList<Integer> cumusCounts = new ArrayList<>();
        ArrayList<Double> cumusProbs = new ArrayList<>();
        ArrayList<Double> cumusPoints = new ArrayList<>();

        StatisticsPackage p = new StatisticsPackage();

        //pivot = max - min / 10;
        int bin_size = cumusPoints_for_reference.size();

        for (int b = 0; b < bin_size; b++) {
            int count = 0;
            double threshold = cumusPoints_for_reference.get(b);

            cumusPoints.add(threshold);
            for (int i = 0; i < sample.size(); i++) {
                //if (sample.get(i) <= (threshold)) - we can lost precision here.
                if (threshold - sample.get(i) > -0.0001) {
                    count++;
                }
            }
            cumusCounts.add(count);
            cumusProbs.add((double) count / (double) sample.size());
        }
        /*for (int b = 0; b < bin_size; b++) {
            System.out.println(cumusCounts.get(b));
        }*/
        p.cumusCounts = cumusCounts;
        p.cumusPoints = cumusPoints;
        p.cumusProbs = cumusProbs;
        return p;
    }

    Random r = new Random();
    public double GaussianPDF(double x, double mean, double variance) {
        return (Math.exp(-(((x - mean) * (x - mean)) / ((2.0 * variance)))))*(1.0 / (Math.sqrt(variance * 2.0 * Math.PI)));
    }
    public double sliceSampling_stepout_procedure(double x_0, double y, double mean, double variance) {
        
        double w = 1.0;
        double m = 100 ;
        
        double u = r.nextDouble();        
        double L = x_0 - w*u;
        double R = L + w;
        double v = r.nextDouble();
        int J = (int) Math.floor((m*v));
        int K = (int) ((m-1)-J);
        
        //double y = GaussianPDF(x_0, mean, variance);
        while(J>0 && y < GaussianPDF(L, mean, variance)) {
            L = L - w;
            J = J-1;
        }
        while( K>0 && y<GaussianPDF(R, mean, variance)) {
            R = R + w;
            K = K - 1;
        }
        
        //shrinkage procedure
        return shrinkageProcedure(L, R, y, x_0, mean, variance);                
    }
    public double sliceSampling_doubling_procedure(double x_0, double y, double mean, double variance) {
        
        double w = 0.5;
        double p = Math.pow(2.0, 10) ;
        
        double u = r.nextDouble();        
        double L = x_0 - w*u;
        double R = L + w;
        
        int K = (int) p;
        
        while( K>0 && (y<GaussianPDF(L, mean, variance) || y<GaussianPDF(R, mean, variance))) {
            double v = r.nextDouble();
            if(v<0.5) {
                L = L - (R - L);
            } else {
                R = R + (R - L);
            }
            K--;
        }
        
        //shrinkage procedure
        return testProcedure(L, R, y, x_0, mean, variance);                
    }
    
    public double testProcedure(double L, double R, double y, double x_0, double mean, double variance) {
        while(true) {
            double u = r.nextDouble();
            double x_1 = L + u*(R-L);
            if(y < GaussianPDF(x_1,mean, variance)) {
                return x_1;
            }
            if(x_1<x_0) {
                L = x_1;
            }
            else {
                R = x_1;
            }
            //System.out.println(x_1);
        }
    }
    
    
    public double shrinkageProcedure(double L, double R, double y, double x_0, double mean, double variance) {
        while(true) {
            double u = r.nextDouble();
            double x_1 = L + u*(R-L);
            if(y < GaussianPDF(x_1,mean, variance)) {
                return x_1;
            }
            if(x_1<x_0) {
                L = x_1;
            }
            else {
                R = x_1;
            }
            //System.out.println(x_1);
        }
    }
    
    public double SampleVariance(ArrayList<Double> list, double mean) {
        double sum = 0;
        for(int i = 0; i < list.size(); i++) {
            sum+=Math.pow(list.get(i)-mean, 2.0);
        }
        sum = sum/((double) list.size());
        return sum;
    }
    
    public double SampleMean(ArrayList<Double> list) {
        double sum = 0;
        for(int i = 0; i < list.size(); i++) {
            sum+=list.get(i);
        }
        return sum/((double) list.size());
    }
    
    public static void main(String args[]) {
        double mean = 0;
        double variance = 5;
        SliceSampler sampler = new SliceSampler();
        double u = 0.0;
        ArrayList<Double> list = new ArrayList<>();
        double rangeMin = 0.0;
        for(int i = 0; i < 10000; i++) {
            double rangeMax = sampler.GaussianPDF(u, mean, variance);
            double randomValue = rangeMin + (rangeMax - rangeMin) * sampler.r.nextDouble();
            
            u = sampler.sliceSampling_stepout_procedure(u, randomValue, mean, variance);
            list.add(u);
        }
        for(int i = 0; i < 100; i++) {
            //System.out.println(list.get(i));
        }
        System.out.println(sampler.SampleMean(list));
        System.out.println(sampler.SampleVariance(list, sampler.SampleMean(list)));
        
        StatisticsPackage p_empirical = sampler.cumuculativedistribution(list);
        
        if (1 == 1) {
            Random r = new Random();
            ArrayList<Double> another_list_of_means = new ArrayList<>();
            for (int i = 0; i < list.size(); i++) {
                double number = r.nextGaussian() * (Math.sqrt(variance)) + mean;
                another_list_of_means.add(number);
            }
            StatisticsPackage p_reference = sampler.cumuculativedistribution(another_list_of_means);

            for (int i = 0; i < p_empirical.cumusPoints.size(); i++) {
                System.out.println(p_empirical.cumusPoints.get(i)+"~"+p_reference.cumusPoints.get(i)+"~"+p_empirical.cumusProbs.get(i) + "~" + p_reference.cumusProbs.get(i));
            }
        }
        
    }
}
