
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Queue;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.Map.Entry;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import com.csvreader.CsvWriter;

public class MainClass {

	//public static HashMap<Integer, ArrayList<Double>> gene = new HashMap<Integer, ArrayList<Double>>();
	public static HashMap<Integer,Integer> groundtruth= new HashMap<Integer,Integer>();

	static HashMap <Integer,ArrayList<Double> > all_datapoints  = new HashMap<>();
	static HashMap <Integer, ArrayList<Double>> centroids  = null;
	static HashMap <Integer, ArrayList<Double>> prev_centroids  = null;
	static HashMap<Integer, ArrayList<ArrayList<Double>>> resultMap = null;
	public static HashMap<Integer,Integer> last_Iter_cluster = new HashMap<Integer, Integer>();

	public static ArrayList<ArrayList<Integer>> dbScanCluster= new ArrayList<ArrayList<Integer>>();
	public static ArrayList<ArrayList<Integer>> clusterList= new ArrayList<ArrayList<Integer>>();
	public static ArrayList<Integer> visitedList =new ArrayList<Integer>();
	public static Map<Integer,Integer> gene_cluster_dbscan= new HashMap<Integer, Integer>();

	public static int Hrclevel=0;
	public static double DBScanEps=0;
	public static int DbScanMinpts=0;
	public static String fileName=null; 

	static int k_centroids = 5;
	//public static double[][] distanceMatrix;


	public static void initialize()
	{
		centroids  = new HashMap<>();


		File configFile = new File("config.properties.txt");
		BufferedReader reader;
		try {


			FileReader reader1 = new FileReader(configFile);
			Properties props = new Properties();
			props.load(reader1);

			k_centroids=Integer.parseInt(props.getProperty("No_of_centroids"));

			
			Hrclevel=Integer.parseInt(props.getProperty("Level"));
			DBScanEps= Double.parseDouble(props.getProperty("eps"));
			DbScanMinpts=Integer.parseInt(props.getProperty("MinPts"));
			fileName=props.getProperty("fileName");
			
			
			
			reader = new BufferedReader(new FileReader(fileName));

			int a =1;
			String line = null;
			while ((line = reader.readLine()) != null) 
			{	
				String str_arr[] = line.split("\t");


				ArrayList<Double> exp= new ArrayList<Double>();

				for(int i = 2; i<str_arr.length; i++)
				{

					exp.add(Double.parseDouble(str_arr[i]));
				}

				all_datapoints.put(Integer.parseInt(str_arr[0]), exp);

				groundtruth.put(Integer.parseInt(str_arr[0]),Integer.parseInt(str_arr[1]));

			}

			HashMap <Integer,ArrayList<Double> > all_datapoints1  = new HashMap<>();
			double[] temp=null;
			int size=all_datapoints.size();
			for (Integer a1 : all_datapoints.keySet()) {
				temp= new double[all_datapoints.get(a1).size()];
				break;
			}

			for (Integer a1 : all_datapoints.keySet()) {
				ArrayList<Double> t= new ArrayList<>();
				t=(ArrayList<Double>) all_datapoints.get(a1).clone();
				for(int i=0;i<t.size();i++)
				{
					temp[i]+=t.get(i);
				}
			}
			for (Integer a1 : all_datapoints.keySet()) {
				
				ArrayList<Double> t= new ArrayList<>();
				t=(ArrayList<Double>) all_datapoints.get(a1).clone();
				ArrayList<Double> t1= new ArrayList<>();
				for(int i=0;i<t.size();i++)
				{
					   t1.add(t.get(i)-(temp[i]/(double)size));
				}
				all_datapoints1.put(a1, t1);
			}

			all_datapoints.clear();
			all_datapoints=(HashMap<Integer, ArrayList<Double>>) all_datapoints1.clone();

			for(int i=1; i<=k_centroids; i++)
			{
				String centroid = props.getProperty(String.valueOf(i));
				centroids.put(i,(ArrayList<Double>) all_datapoints.get(Integer.parseInt(centroid)).clone());
			}
			
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}



	}


	public static double findEuclidianDistance(ArrayList<Double> exp, ArrayList<Double> centroid)
	{
		double Dist = 0;
		for(int i=0; i< centroid.size();i++)
		{
			double val= exp.get(i)-centroid.get(i);
			Dist+= (val*val); 
		}
		return Math.sqrt(Dist);


	}



	//-----------------KMeans clustring---------------------------------------------
	public static void Kmeans()
	{


		boolean flag = true;

		int loops =0;

		int count =0;
		while(flag == true)
			//while(count<4)
		{
			Kmeans1(all_datapoints,centroids);
			flag = calcNewCentroids(resultMap);
			System.out.println();
			System.out.println("########################################################  loop count "+ ++loops);
			count++;
		}


		//		System.out.println("--------------------------Final Clusters------------------------------");
		//		for(Entry<Integer, Integer> entry : last_Iter_cluster.entrySet())
		//		{
		//			System.out.println("gene_id  "+entry.getKey()+"   cluster_no  "+entry.getValue());
		//		}

		//	StringBuilder sb= new StringBuilder();
		//		
		//		
		//		for(int i=1;i<=all_datapoints.size();i++)
		//		{
		//			sb.append(i+"\t");
		//			if(last_Iter_cluster.get(i)!=null)
		//			{
		//				sb.append(last_Iter_cluster.get(i)+"\t");
		//			}
		//			else
		//			{
		//				sb.append(-1+"\t");
		//			}
		//			ArrayList<Double> arr= new ArrayList<Double>();
		//			arr= (ArrayList<Double>) all_datapoints.get(i).clone();
		//			for (Double double1 : arr) {
		//				sb.append(double1+"\t");
		//			}
		//			
		//			sb.trimToSize();
		//			sb.append("\n");
		//		}
		//		try {
		//		File file = new File("kmeans.txt");
		//		if (!file.exists()) {
		//			
		//				file.createNewFile();
		//			
		//		}
		//		FileWriter fw = new FileWriter(file.getAbsoluteFile());
		//		BufferedWriter bw = new BufferedWriter(fw);
		//		bw.write(sb.toString());
		//		bw.close();
		//		} catch (IOException e) {
		//			// TODO Auto-generated catch block
		//			e.printStackTrace();
		//		}


		String outputFile = "kmeans.csv";

		try {
			// use FileWriter constructor that specifies open for appending
			CsvWriter csvOutput = new CsvWriter(new FileWriter(outputFile, true), ',');



			for(int i=1;i<=all_datapoints.size();i++)
			{
				csvOutput.write(String.valueOf(i));
				if(last_Iter_cluster.get(i)!=null)
				{

					csvOutput.write(String.valueOf(last_Iter_cluster.get(i)));
				}
				else
				{
					csvOutput.write("-1");
				}
				ArrayList<Double> arr= new ArrayList<Double>();
				arr= (ArrayList<Double>) all_datapoints.get(i).clone();
				for (Double double1 : arr) {
					//sb.append(double1+"\t");
					csvOutput.write(String.valueOf(double1));
				}
				csvOutput.endRecord();

			}



			csvOutput.close();
		} catch (IOException e) {
			e.printStackTrace();
		}


		System.out.println();
		double jaccard= jaccardCo(groundtruth, last_Iter_cluster);
		System.out.println("jaccard  :"+jaccard);
		double coorelation =coRelation(last_Iter_cluster, all_datapoints);
		System.out.println("CoOrelation  :"+coorelation);

	}

	public static void Kmeans1(HashMap<Integer, ArrayList<Double>> geneData, HashMap<Integer, ArrayList<Double>> centroids)
	{

		resultMap= new HashMap<>();
		last_Iter_cluster = new HashMap<Integer, Integer>();

		for (int gene:geneData.keySet()) 
		{
			ArrayList<Double> expression= new ArrayList<Double>();
			expression= geneData.get(gene);
			double currentMin=Double.POSITIVE_INFINITY;
			int clusterNo=0;
			for(int centro:centroids.keySet())
			{
				ArrayList<Double> ctr= new ArrayList<Double>();
				ctr= centroids.get(centro);
				double mindiff= findEuclidianDistance(expression, ctr);
				if(mindiff<currentMin)
				{
					currentMin= mindiff;
					clusterNo=centro;
				}
			}

			if(resultMap.containsKey(clusterNo))
			{
				ArrayList<ArrayList<Double>> t= new ArrayList<>();
				t = resultMap.get(clusterNo);
				t.add(geneData.get(gene));
				resultMap.put(clusterNo, t);

				last_Iter_cluster.put(gene, clusterNo);
			}
			else
			{

				ArrayList<ArrayList<Double>> t= new ArrayList<>();
				t.add(geneData.get(gene));
				resultMap.put(clusterNo, t);

				last_Iter_cluster.put(gene, clusterNo);
			}

		}


		//		int abc =0;
		//		for(Entry<Integer, ArrayList<ArrayList<Double>>> entry : resultMap.entrySet())
		//		{
		//
		//			System.out.println();
		//			System.out.println("Cluster ------  "+entry.getKey());
		//			ArrayList<ArrayList<Double>>  temp = entry.getValue();
		//
		//			for (ArrayList<Double> arrayList : temp) 
		//			{
		//				System.out.print(++abc+"     ");
		//				for (Double  d : arrayList) 
		//				{
		//					System.out.print("   "+d);
		//				}
		//				System.out.println();
		//			}
		//
		//		}



	}

	public static boolean calcNewCentroids(HashMap<Integer, ArrayList<ArrayList<Double>>> resultMap2)
	{
		boolean flag = true;
		prev_centroids = new HashMap<Integer, ArrayList<Double>>();
		prev_centroids = (HashMap<Integer, ArrayList<Double>>) centroids.clone();
		centroids  = new HashMap<>();
		int len =0;

		//		System.out.println(" Previous Centroids ---------------------------");
		//
		//		for(Entry<Integer, ArrayList<Double>> entry : prev_centroids.entrySet())
		//		{
		//			System.out.print("The Key is  "+entry.getKey()+"   ");
		//
		//			for (Double arrayList : entry.getValue()) 
		//			{
		//				System.out.print(arrayList+"  ");
		//			}
		//			System.out.println();
		//		}


		for(Entry<Integer, ArrayList<ArrayList<Double>>> entry : resultMap2.entrySet())
		{

			ArrayList<ArrayList<Double>> list = entry.getValue();


			//ArrayList<Integer> temp = new ArrayList<Integer>();
			len = entry.getValue().get(0).size();

			int l = entry.getValue().size();

			double temp[] =  new double[len];

			for (ArrayList<Double> arrayList : list) {

				Object[] first = arrayList.toArray();
				for(int j =0 ; j< first.length; j++ )
				{
					//double d = (double) first[j];
					temp[j]=temp[j]+ (double)first[j];
				}


			}


			ArrayList<Double> a = new ArrayList<>();	

			for(int k =0 ; k< len; k++ )
			{
				a.add(temp[k]/(double)l);
			}



			centroids.put(entry.getKey(), a);

		}

		int c =0;



		int no_of_centroid=prev_centroids.size();
		//		System.out.println("Number of Previous Clusters   " +no_of_centroid);
		//		System.out.println("Number of Current Clusters   " +centroids.size());

		//System.out.println("**********************************************");

		for(Entry<Integer, ArrayList<Double>> e : prev_centroids.entrySet())
		{
			//			System.out.println();
			//			System.out.println(e.getKey()+"  ^^^^^^^^^^^^^^^^^^^^");

			int key = e.getKey();
			ArrayList<Double> prev = new ArrayList<Double>();
			prev = e.getValue();

			//			System.out.println("Previous");
			//			for (Double double1 : prev) {
			//				System.out.print(double1.doubleValue()+"--");
			//			}
			//			System.out.println();

			if(centroids.containsKey(key))
			{
				ArrayList<Double> curr = new ArrayList<Double>();
				curr = centroids.get(key);


				System.out.println("Current");
				for (Double double2 : curr) {
					System.out.print(double2.doubleValue()+"\t");
				}



				Object[] c_a = curr.toArray();;
				Object[] p_a = prev.toArray();

				if(Arrays.deepEquals(c_a,  p_a))
				{
					c++;
				}

			}


		}

		if(c == no_of_centroid) 
		{
			return false ;
		}
		else
		{
			return true;	
		}

	}


	//-----------------------------------------------------------------------------



	//-----------Heirarchial clustring-------------------------------------
	public static void HRC(HashMap<Integer, ArrayList<Double>> geneData)
	{
		ArrayList<ArrayList<String>> main_cluster = new ArrayList<ArrayList<String>>();
		ArrayList<String> cluster = new ArrayList<String>();
		HashMap<String, Double> distanceMap= new HashMap<>();
		for(int geneid:geneData.keySet())
		{
			cluster.add(String.valueOf(geneid));
		}
//		System.out.println("--------------------------------Level 1----------------------------------"+cluster.size());
//		for (String string : cluster) {
//			System.out.println(string);
//		}

		main_cluster.add(cluster);

		int size= geneData.size();
		ArrayList<Double> tempExp1= new ArrayList<Double>();
		ArrayList<Double> tempExp2= new ArrayList<Double>();

		double min= Double.POSITIVE_INFINITY;
		String KeyMin= null;
		String cluster1=null;
		String cluster2= null;
		for(int i=1;i<=size;i++)
		{
			tempExp1= geneData.get(i);

			for(int j=i+1;j<=size;j++)
			{

				tempExp2= geneData.get(j);
				double dist= findEuclidianDistance(tempExp1, tempExp2);
				String key= String.valueOf(i)+"_"+String.valueOf(j);
				distanceMap.put(key, dist);
				String key1= String.valueOf(j)+"_"+String.valueOf(i);
				distanceMap.put(key1, dist);
				if(dist<min)
				{
					min=dist;
					KeyMin=key;
					cluster1= String.valueOf(i);
					cluster2= String.valueOf(j);

				}

			}
		}

		ArrayList<String> tempCluster= new ArrayList<>();
		for (String string : cluster) {
			if((!string.equals(cluster1))&&(!string.equals(cluster2)))
				tempCluster.add(string);
		}
		tempCluster.add(KeyMin);
//		System.out.println("--------------------------------Level 2----------------------------------"+tempCluster.size());
//		for (String string : tempCluster) {
//			System.out.println(string);
//		}

		main_cluster.add(tempCluster);

		int level=3;
		while(tempCluster.size()>1)
		{
			cluster= new ArrayList<>();
			double iterativeMin= Double.POSITIVE_INFINITY;
			String c1= null;
			String c2= null;
			String c3= null;
			for (int i=0;i<tempCluster.size();i++) {
				for(int j=i+1;j<tempCluster.size();j++)
				{
					String[] s1= tempCluster.get(i).split("_");
					String[] s2= tempCluster.get(j).split("_");
					double minIntercluster=Double.POSITIVE_INFINITY;
					for (String string : s1) {
						for (String string1 : s2) {
							if(minIntercluster>distanceMap.get(string+"_"+string1))
							{
								minIntercluster= distanceMap.get(string+"_"+string1);
							}

						}

					}
					if(iterativeMin>minIntercluster)
					{
						iterativeMin=minIntercluster;
						c1=tempCluster.get(i);
						c2=tempCluster.get(j);
						c3= c1+"_"+c2;

					}



				}

			}

			for (String string : tempCluster) {
				if((!string.equals(c1))&&(!string.equals(c2)))
					cluster.add(string);
			}
			cluster.add(c3);
			main_cluster.add(cluster);

			tempCluster.clear();
			tempCluster=(ArrayList<String>) cluster.clone();
//			System.out.println("--------------------------------Level "+level+++"----------------------------------"+tempCluster.size());
//			for (String string : tempCluster) {
//				System.out.println(string);
//			}

		}

		last_Iter_cluster= new HashMap<Integer, Integer>();
		ArrayList<String> show= new ArrayList<String>();
		show= main_cluster.get(geneData.size()-Hrclevel);
		int index=1;
		for (String string : show) {
			String[] spl= string.split("_");
			for (String string2 : spl) {
				last_Iter_cluster.put(Integer.parseInt(string2), index);
			}
			index++;
		}

		String outputFile = "HRC.csv";

		try {
			// use FileWriter constructor that specifies open for appending
			CsvWriter csvOutput = new CsvWriter(new FileWriter(outputFile, true), ',');



			for(int i=1;i<=all_datapoints.size();i++)
			{
				csvOutput.write(String.valueOf(i));
				if(last_Iter_cluster.get(i)!=null)
				{

					csvOutput.write(String.valueOf(last_Iter_cluster.get(i)));
				}
				else
				{
					csvOutput.write("-1");
				}
				ArrayList<Double> arr= new ArrayList<Double>();
				arr= (ArrayList<Double>) all_datapoints.get(i).clone();
				for (Double double1 : arr) {
					//sb.append(double1+"\t");
					csvOutput.write(String.valueOf(double1));
				}
				csvOutput.endRecord();

			}



			csvOutput.close();
		} catch (IOException e) {
			e.printStackTrace();
		}


		double jaccard= jaccardCo(groundtruth, last_Iter_cluster);
		System.out.println("jaccard  :"+jaccard);
		double coorelation =coRelation(last_Iter_cluster, all_datapoints);
		System.out.println("CoOrelation  :"+coorelation);
	}
	//--------------------------------------------------------------------




	//------------Density based cluster--------------------------
	public static void DBSCAN(HashMap<Integer, ArrayList<Double>> geneData)
	{
		ArrayList<Integer>Cluster;
		double eps=DBScanEps;
		int minPts= DbScanMinpts;
		HashMap<Integer,Boolean>Visited=new HashMap<Integer,Boolean>();
		ArrayList<Integer> outlierGenes= new ArrayList<Integer>();
		for(int i:geneData.keySet())
		{
			Visited.put(i, false);
		}

		for(int i: geneData.keySet())
		{
			if(Visited.get(i)==false)
			{
				Visited.put(i, true);
				ArrayList<Integer> Neighbours= new ArrayList<Integer>();
				Neighbours=regionQuery(i,eps,geneData);
				if(Neighbours.size()>=minPts)
				{
					Cluster= new ArrayList<Integer>();
					ArrayList<Integer> expndCluster= expandCluster(i,Neighbours,Cluster,eps,minPts,geneData, Visited);
					dbScanCluster.add(expndCluster);
				}
				else
				{
					outlierGenes.add(i);
				}
			}
		}
		System.out.println();
		last_Iter_cluster= new HashMap<Integer, Integer>();
		int index=1;

		for (ArrayList<Integer> ar1 : dbScanCluster) {
			for (Integer integer : ar1) {
				last_Iter_cluster.put(integer, index);

			}
			index++;

		}

		String outputFile = "dbscan.csv";

		try {
			// use FileWriter constructor that specifies open for appending
			CsvWriter csvOutput = new CsvWriter(new FileWriter(outputFile, true), ',');



			for(int i=1;i<=all_datapoints.size();i++)
			{
				csvOutput.write(String.valueOf(i));
				if(last_Iter_cluster.get(i)!=null)
				{

					csvOutput.write(String.valueOf(last_Iter_cluster.get(i)));
				}
				else
				{
					csvOutput.write("-1");
				}
				ArrayList<Double> arr= new ArrayList<Double>();
				arr= (ArrayList<Double>) all_datapoints.get(i).clone();
				for (Double double1 : arr) {
					//sb.append(double1+"\t");
					csvOutput.write(String.valueOf(double1));
				}
				csvOutput.endRecord();

			}



			csvOutput.close();
		} catch (IOException e) {
			e.printStackTrace();
		}


		double jaccard= jaccardCo(groundtruth, last_Iter_cluster);
		System.out.println("jaccard  :"+jaccard);
		double coorelation =coRelation(last_Iter_cluster, all_datapoints);
		System.out.println("CoOrelation  :"+coorelation);

	}

	public static ArrayList<Integer> regionQuery(int geneId, double eps,HashMap<Integer, ArrayList<Double>> geneData)
	{
		ArrayList<Integer> neighbour=new ArrayList<Integer>();
		ArrayList<Double> exp1=geneData.get(geneId);

		for(int geneId1:geneData.keySet())
		{
			if(geneId!=geneId1)
			{
				ArrayList<Double> exp2=geneData.get(geneId1);
				Double dist=findEuclidianDistance(exp1, exp2);
				if(dist<=eps)
				{
					neighbour.add(geneId1);
				}
			}
		}
		return neighbour;
	}

	public static ArrayList<Integer> expandCluster(int geneId,ArrayList<Integer> neighbour,ArrayList<Integer> cluster,double eps,
			int minPts,HashMap<Integer, ArrayList<Double>> geneData,HashMap<Integer,Boolean>Visited )
			{

		HashSet<Integer> s=new HashSet<Integer>();
		s.addAll(neighbour);
		Queue<Integer> seeds = new LinkedList<Integer>();
		seeds.addAll(neighbour);
		cluster.add(geneId);

		while(!seeds.isEmpty())
		{
			int nbrgene=(Integer) seeds.poll();
			if(Visited.get(nbrgene)==false)
			{
				Visited.put(nbrgene, true);
				ArrayList<Integer> nbrNeighbour=new ArrayList<Integer>();
				nbrNeighbour=regionQuery(nbrgene, eps, geneData);
				if(nbrNeighbour.size()>=minPts)
				{
					seeds.addAll(nbrNeighbour);
					//s.addAll(nbrNeighbour);
				}
			}
			if(dbScanCluster.size()>0)
			{
				for(ArrayList<Integer> arr:dbScanCluster)
				{
					if(!arr.contains(nbrgene))
					{
						if(!cluster.contains(nbrgene))
						{
							cluster.add(nbrgene);
						}
					}
				}
			}
			else
			{
				if(!cluster.contains(nbrgene))
				{
					cluster.add(nbrgene);
				}
			}

			//			int countlist=0;
			//			for(ArrayList<Integer> list:dbScanCluster)
			//			{
			//				if(!list.contains(nbrgene))
			//				{
			//				countlist++;	
			//				}
			//			}
			//			if(countlist==dbScanCluster.size())
			//			{
			//				if(!cluster.contains(nbrgene))
			//					cluster.add(nbrgene);
			//			}
		}
		return cluster;

			}
	//---------------------------------------------------------------


	//--------------------External Index----------------------------------------

	public static double jaccardCo(HashMap<Integer,Integer> groundtruth,HashMap<Integer,Integer> ourcluster)
	{
		double jaccard= 0;
		double M11 = 0;		
		double M01 = 0;     
		double M10 = 0;     
		int size = ourcluster.size();
		int[][] ourClusterMat = new int[size][size];
		int[][] grnClusterMat = new int[size][size];

		for(int i = 0; i < size; i++) {
			for(int j = i; j < size; j++) {
				if(ourcluster.get(i)==null||ourcluster.get(j)==null)
				{
					ourClusterMat[i][j] = ourClusterMat[j][i] = 0;
				}
				else
				{
					if(ourcluster.get(i) == ourcluster.get(j))
						ourClusterMat[i][j] = ourClusterMat[j][i] = 1;
					else
						ourClusterMat[i][j] = ourClusterMat[j][i] = 0;
				}

				if(groundtruth.get(i) == groundtruth.get(j))
					grnClusterMat[i][j] = grnClusterMat[j][i] = 1;
				else
					grnClusterMat[i][j] = grnClusterMat[j][i] = 0;
			}
		}


		for(int i = 0; i < size; i++) {
			for(int j = 0; j < size; j++) {
				if(ourClusterMat[i][j] == grnClusterMat[i][j]) {
					if(ourClusterMat[i][j] == 1)
						M11++;
				}
				else {
					if(ourClusterMat[i][j] == 1 && grnClusterMat[i][j] == 0)
						M10++;
					else if(ourClusterMat[i][j] == 0 && grnClusterMat[i][j] == 1)
						M01++;
				}
			}
		}
		jaccard= M11/(M11+M10+M01);
		return jaccard;
	}

	//------------------Internal index-----------------------------------------

	public static double coRelation(Map<Integer,Integer> gene_cluster_dbscan, HashMap <Integer,ArrayList<Double>> gene)
	{
		double coRel =0.0;
		int size = gene.size();
		double [] D = getDistMat(gene);
		double [] C = new double [size*size];
		int cnt =0;
		for(int i=1;i<=size;i++)
			for(int j=1;j<=size;j++)
			{
				if(gene_cluster_dbscan.get(i)==null||gene_cluster_dbscan.get(j)==null)
				{
					C[cnt] =0;
				}
				else
				{
					if(gene_cluster_dbscan.get(i) == gene_cluster_dbscan.get(j))
						C[cnt] = 1;
					else
						C[cnt] =0;
				}
				cnt++;

			}
		coRel= new PearsonsCorrelation().correlation(D, C);
		//return Math.abs(coRel * ((size*(size - 1)/2)));
		return coRel;
	}

	public static double[] getDistMat(HashMap <Integer,ArrayList<Double>> gene)
	{
		int size = gene.size();
		double[] D_mat= new double[size*size];
		int count = 0;

		ArrayList<Double> tempExp1= new ArrayList<Double>();
		ArrayList<Double> tempExp2= new ArrayList<Double>();

		for(int i=1;i<=size;i++)
		{
			tempExp1= gene.get(i);
			double min= Double.POSITIVE_INFINITY;
			for(int j=1;j<=size;j++)
			{
				if(i==j)
				{
					D_mat[count] = 0;
				}
				else
				{

					tempExp2= gene.get(j);
					D_mat[count]= findEuclidianDistance(tempExp1, tempExp2);


				}
				count++;

			}
		}
		return D_mat;
	}

	//--------------------Main-------------------------------------------	
	public static void main(String[] args) {


		initialize();
		Kmeans();
		HRC(all_datapoints);
		DBSCAN(all_datapoints);
	
	}





}
