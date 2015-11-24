package edu.asu.irs13;

import java.io.*;
import java.util.*;

import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.Term;
import org.apache.lucene.index.TermDocs;
import org.apache.lucene.index.TermEnum;
import org.apache.lucene.store.FSDirectory;

public class IR_PrjPrt2_PageRank_by_Priyank {
	
	/*----------------------------------------------------------------------Code from the LinkAnalysis.java file------------------------------------*/
	public static final String linksFile = "IntLinks.txt";
	public static final String citationsFile = "IntCitations.txt";
	public static int numDocs = 25053;

	private int[][] links;
	private int[][] citations;
	public static double resetValue = 1.0/25054.0;
	
	public IR_PrjPrt2_PageRank_by_Priyank()
	{
		try
		{
			// Read in the links file
			links = new int[numDocs][];
			BufferedReader br = new BufferedReader(new FileReader(linksFile));
			String s = "";
			while ((s = br.readLine())!=null)
			{
				String[] words = s.split("->"); // split the src->dest1,dest2,dest3 string
				int src = Integer.parseInt(words[0]);
				if (words.length > 1 && words[1].length() > 0)
				{
					String[] dest = words[1].split(",");
					links[src] = new int[dest.length];
					for (int i=0; i<dest.length; i++)
					{
						links[src][i] = Integer.parseInt(dest[i]);
					}
				}
				else
				{
					links[src] = new int[0];
				}
			}
			br.close();
			
			// Read in the citations file
			citations = new int[numDocs][];
			br = new BufferedReader(new FileReader(citationsFile));
			s = "";
			while ((s = br.readLine())!=null)
			{
				String[] words = s.split("->"); // split the src->dest1,dest2,dest3 string
				int src = Integer.parseInt(words[0]);
				if (words.length > 1 && words[1].length() > 0)
				{
					String[] dest = words[1].split(",");
					citations[src] = new int[dest.length];
					for (int i=0; i<dest.length; i++)
					{
						citations[src][i] = Integer.parseInt(dest[i]);
					}
				}
				else
				{
					citations[src] = new int[0];
				}

			}
			br.close();
		}
		catch(NumberFormatException e)
		{
			System.err.println("links file is corrupt: ");
			e.printStackTrace();			
		}
		catch(IOException e)
		{
			System.err.println("Failed to open links file: ");
			e.printStackTrace();
		}
	}
	
	public int[] getLinks(int docNumber)
	{
		return links[docNumber];
	}
	
	public int[] getCitations(int docNumber)
	{
		return citations[docNumber];
	
	}
	/*-----------------------------------------------------------Code from the LinkAnalysis.java file -----------------------------------------------------------------------*/
	
	//Start of the Main function
	public static void main(String[] args) throws IOException {
		
		//Initializing the numDocs value and creating a new variable l to access the links and citations retrieved above.
		IR_PrjPrt2_PageRank_by_Priyank.numDocs = 25054;
		IR_PrjPrt2_PageRank_by_Priyank l = new IR_PrjPrt2_PageRank_by_Priyank();
		
		//Initializing the value of 'w'
		double w = 0.4;
		
		//Declaring the rank_current,Previous Rank and M array.
		double[] rank_curr = new double[IR_PrjPrt2_PageRank_by_Priyank.numDocs];
		double[] M = new double[IR_PrjPrt2_PageRank_by_Priyank.numDocs];
		double[] rank_prev = new double[IR_PrjPrt2_PageRank_by_Priyank.numDocs];
		
		//initializing the LargestError and the Threshold value.
		double LargestError = 1;
		double Threshold = 0.0000001;
		
		//Initializing the array with initial values
		rank_curr = CreateInitialMatrix(IR_PrjPrt2_PageRank_by_Priyank.numDocs);
		
		/*---------------------------------------------Convergence loop to calculate the global PageRank after going through the power iterations----------------------------*/ 
		
		do{
			rank_prev = equalToArray(rank_curr);
			
			for(int docs = 0; docs < IR_PrjPrt2_PageRank_by_Priyank.numDocs; docs++){
				M = compute(docs,l);
				rank_curr[docs] = MatrixMultiply(M,rank_curr);
			}
			
			rank_curr = NormalizeMatrix(rank_curr);
			LargestError = findingLargest(rank_prev,rank_curr);
			
			
		}while(LargestError > Threshold);
		
		/*----------------------------------------End of Convergence loop------------------------------------------------------------*/
		
		
		/*-------------------------------------MIN-MAX Normalization---------------------------*/
		
		double min = 0;
		double max = 0.00000000000000000000001;
		for(int row = 0; row < IR_PrjPrt2_PageRank_by_Priyank.numDocs; row++){
			if(rank_curr[row] < min)
				min = rank_curr[row];
			if(rank_curr[row] > max)
				max = rank_curr[row];
		}
		
		double difference = max - min;
		
		for(int row = 0; row < IR_PrjPrt2_PageRank_by_Priyank.numDocs; row++)
			rank_curr[row] = (rank_curr[row] - min)/difference;

	
		/*-----------------------------END of Min-Max--------------------------------------------------------*/
		
		/*---------------------------PAGE RANK VARIABLE--------------------------------------*/
		//HashMap<Integer, Double> finalPageRank = new HashMap<Integer, Double>();
		
		/*-----------------------------------------------------------Start of Code from Project Part 1-----------------------------------------------------------------------------*/		
	
		int countDocofTerm = 0;
		double TempIDFValue = 0;
		
		//HashMap to Store the IDF per term.
		HashMap<String,Double> TermIDF = new HashMap<String, Double>();
		
		//Variable to store the all the terms being indexed.
		HashMap<String,Boolean> TermData = new HashMap<String,Boolean>();
		
		//Variable for constructing the Document Term Matrix.
		//It is a HashMap of Integer(which here shall indicate the document number[staring from 0] as the Key
		//AND
		//another HashMap as the Value. This inner HashMap would consist of Term and Frequency as the Key and Value pair.
		HashMap<Integer,HashMap<String,Integer>> dtMatrix = new HashMap<Integer,HashMap<String,Integer>>();

		//Similar Variable for IDF
		HashMap<Integer,HashMap<String,Double>> dtMatrix_TF_IDF = new HashMap<Integer,HashMap<String,Double>>();
		
		//Variable for storing the Magnitude of the Document Vectors during Task 2,i.e. TF/IDF. Key being the Document ID and Double is the Value.
		HashMap<Integer,Double> MagDocVec_TF_IDF = new HashMap<Integer,Double>();
		
		
	

		// the IndexReader object is the main handle that will give you 
		// all the documents, terms and inverted index
		/*---Copied from the SeachFiles.Java----*/
		IndexReader r = IndexReader.open(FSDirectory.open(new File("index")));

		//get Total number of Documents
		int totDocs = r.maxDoc();


		TermEnum t = r.terms();
		t = r.terms();

		//Loop to go through every term and construct it individual Document-Term Matrix
		while(t.next()){
			if(!(TermData.containsKey(t.term().text())))
				TermData.put(t.term().text(), true);

			Term te = new Term("contents", t.term().text());
			TermDocs td = r.termDocs(te);

			while(td.next()){
				HashMap<String,Integer>	docRow = new HashMap<String,Integer>();

				if(dtMatrix.containsKey(td.doc())){
					docRow = dtMatrix.get(td.doc());
				}

				docRow.put(t.term().text(), td.freq());
				dtMatrix.put(td.doc(), docRow);
			}
		}

		
		t = r.terms();
		
		//Loop to go through every term and construct it individual Document-Term Matrix
		while(t.next())
		{
			countDocofTerm = 0;

			Term te1 = new Term("contents", t.term().text());
			TermDocs td1 = r.termDocs(te1);

			while(td1.next()){
				countDocofTerm++;
			}

			//Adding the IDF to the HashMap against the Term.
			TermIDF.put(t.term().text(), Math.log((double) totDocs/((double) countDocofTerm)));
		}





		Set outer = dtMatrix.entrySet();
		Iterator loop = outer.iterator();

		while(loop.hasNext()) {
			Map.Entry outerMap = (Map.Entry)loop.next();

			HashMap<String,Integer>	docRow = (HashMap<String,Integer>) outerMap.getValue(); 
			HashMap<String,Double>	docRow_TF_IDF = new HashMap<String,Double>();

			Set inner = docRow.entrySet();
			Iterator innerLoop = inner.iterator();

			while(innerLoop.hasNext()) {
				Map.Entry<String,Integer> innerMap = (Map.Entry<String,Integer>)innerLoop.next();
				TempIDFValue =(double) innerMap.getValue()*TermIDF.get((String) innerMap.getKey());
				docRow_TF_IDF.put((String) innerMap.getKey(), TempIDFValue);

				if(MagDocVec_TF_IDF.containsKey((int)outerMap.getKey()))
					MagDocVec_TF_IDF.put((int)outerMap.getKey(), MagDocVec_TF_IDF.get((int)outerMap.getKey())+Math.pow(TempIDFValue,2));
				else
					MagDocVec_TF_IDF.put((int)outerMap.getKey(), Math.pow(TempIDFValue,2));
			}

			dtMatrix_TF_IDF.put((int) outerMap.getKey(), docRow_TF_IDF);
		}

		

		double SoP_TF_IDF = 0;
		
		
		Scanner sc = new Scanner(System.in);
		String str = "";
		
		System.out.println("Generating the system. Please wait.........................");
		
		System.out.print("query> ");
		while(!(str = sc.nextLine()).equals("")){
			
			String[] terms = str.split("\\s+");
			
			HashMap<String,Double>	docRow_TF_IDF = new HashMap<String,Double>();
			HashMap<Integer,Double> MagQueryVec = new HashMap<Integer,Double>();
			
			
			//Calculating the Magnitude of the Query Vector
			double disQ = 0;
			for(String word: terms){
				if(TermData.containsKey(word))
					disQ = disQ + 1;
			}
			
			disQ = Math.round(Math.sqrt(disQ)*1000000)/1000000;
			
			//Calculating the final TF into a HashMap
			for(int count = 0;count<totDocs;count++){
				SoP_TF_IDF = 0;
				
				docRow_TF_IDF = dtMatrix_TF_IDF.get(count);
				for(String word : terms)
				{
					if(docRow_TF_IDF.containsKey(word)){
						SoP_TF_IDF = SoP_TF_IDF + docRow_TF_IDF.get(word);
					}
				}
	
				/*---------------------------------------------------------End of the code from Project Part 1-----------------------------------------*/
				
				
				/*------------------------Combining the PageRank and the Vector Similarity for the final values----------------------------------------*/
				if(w == 1.0){
					MagQueryVec.put(count, (rank_curr[count]*w));
				}
					
				
				if(SoP_TF_IDF != 0 && w != 1.0){
					MagQueryVec.put(count,SoP_TF_IDF/(disQ * Math.round(Math.sqrt(MagDocVec_TF_IDF.get(count))*100000)/100000));
					MagQueryVec.put(count, MagQueryVec.get(count)*(1-w) + (rank_curr[count]*w));	//Actual Calculation based on the formula
				}
				
				/*------------------------End of PageRank and Vector combination-----------------------------------------------------------------*/
			}
			
			
			
			// calling a sorting function
			HashMap<Integer,Double> FinalSortedList = HashSort(MagQueryVec);
			
			
			/*-------3. Iterating through the Sorted List and Displaying the Ranked documents------*/
			Set RankedList = FinalSortedList.entrySet();
			loop = RankedList.iterator();
			
			int i = 0;
			
			
			while(loop.hasNext()) {
				Map.Entry finalMap = (Map.Entry)loop.next();
				
				//Displaying the final rank values
				System.out.println("Doc: "+finalMap.getKey()+"]");
				
				if(i == 9)
					break;
				i++;
			}		 
			System.out.print("");
			System.out.print("query> ");
		}
		sc.close();
		
		
		/*------------------------------------------------------------------------------------------*/
		

	}
	
	/*
	 * Initializing the Array function with the reset value
	 */
	public static double[] CreateInitialMatrix(int size){
		double[] iMatrix = new double[size];
		for(int row = 0; row < size; row++)
			iMatrix[row] = 1.0/size;
		return iMatrix;
	}
	
	/*
	 * Function to store value from one array into another
	 */
	public static double[] equalToArray(double[] copyfromArray){
		double[] copyToArray = new double[copyfromArray.length];
		for(int row = 0; row < copyToArray.length; row++)
			copyToArray[row] = copyfromArray[row];
		return copyToArray;
	}
	
	
	/*
	 * The compute Function.
	 * This function computes the M* according to the formula considering the sink nodes computation as well,if any.
	 */
	public static double[] compute(int docNum, IR_PrjPrt2_PageRank_by_Priyank LA){
		double[] result = new double[IR_PrjPrt2_PageRank_by_Priyank.numDocs];
		
		int[] links = LA.getLinks(docNum);
		int[] Cit = LA.getCitations(docNum);
		int[] linksToCount;
		
		//The 'c' variable
		double c = 0.8;
		
		List<Integer> forwardLinks = new ArrayList<Integer>(links.length);
		List<Integer> backwardLinks = new ArrayList<Integer>(links.length);
		
		for (int i=0; i<links.length; i++)
			forwardLinks.add(links[i]);
		for (int i=0; i<Cit.length; i++)
			backwardLinks.add(Cit[i]);
		
		for(int row = 0; row < IR_PrjPrt2_PageRank_by_Priyank.numDocs; row++){
			linksToCount = LA.getLinks(row);
			if(backwardLinks.contains(row)){
				if(linksToCount.length > 1)
					result[row] = (c*(1.0/(double) linksToCount.length)) + ((1.0-c)*resetValue);
				else if(linksToCount.length == 1 && linksToCount[0] == row)
					result[row] = (c) + ((1.0-c)*resetValue);
			}
			else{
				if(linksToCount.length == 0)
					result[row] = c*resetValue + (1.0-c)*resetValue;
				else
					result[row] = (1.0-c)*resetValue;
			}
		}
		
		return result;
		
	}
	
	/*
	 * Function to perform Matrix Multiplication
	 */
	public static double MatrixMultiply(double[] first, double[] second){
		double result = 0;
		for(int row = 0; row < IR_PrjPrt2_PageRank_by_Priyank.numDocs; row++)
				result = result + (first[row]*second[row]);
		
		return result;
	}
	
	/*
	 * Function to Normalize the PageRanks
	 */
	public static double[] NormalizeMatrix(double[] targetVector){
		double L2Norm = 0;
		for(int row = 0; row < IR_PrjPrt2_PageRank_by_Priyank.numDocs; row++)
			L2Norm = L2Norm + targetVector[row];
		for(int row = 0; row < IR_PrjPrt2_PageRank_by_Priyank.numDocs; row++)
			targetVector[row] = targetVector[row]/L2Norm;
		
		return targetVector;
	}
	
	
	/*
	 * Function to find the largest value
	 */
	public static double findingLargest(double[] firstVector, double[] secondVector){
		double largestValue = 0;
		double checkValue = 0;
		for(int row = 0; row < firstVector.length; row++){
			checkValue = Math.abs(firstVector[row] - secondVector[row]);
			if(checkValue >= largestValue)
				largestValue = checkValue;
		}
		
		return largestValue;
	}
	
	
	/*
	 * Function to sort the HashMap w.r.t to the Values.
	 */
	private static HashMap<Integer,Double> HashSort(HashMap<Integer,Double> MGQ) {
		List<Object> TempList = new LinkedList<Object>(MGQ.entrySet());
		
		Collections.sort(TempList, new Comparator<Object>() {
			public int compare(Object Val1, Object Val2) {
				return ((Comparable) ((Map.Entry) (Val2)).getValue()).compareTo(((Map.Entry) (Val1)).getValue());
			}
		});

		//Using the LinkedHashMap, copying the contents into another HashMap while preserving the insertion order
		HashMap<Integer,Double> FinalSortedList = new LinkedHashMap<Integer,Double>();
		
		for (Iterator loop = TempList.iterator(); loop.hasNext();) {
			Map.Entry<Integer,Double> entry = (Map.Entry) loop.next();
			FinalSortedList.put((int)entry.getKey(), (double) entry.getValue());
		}
   
		return FinalSortedList;
	}

}
