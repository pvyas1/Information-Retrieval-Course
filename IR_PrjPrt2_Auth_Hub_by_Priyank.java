package edu.asu.irs13;

import java.io.*;
import java.util.*;
import java.util.concurrent.TimeUnit;

import org.apache.lucene.index.*;
import org.apache.lucene.store.*;


public class IR_PrjPrt2_Auth_Hub_by_Priyank {

/*----------------------------------------------------------------------Code from the LinkAnalysis.java file------------------------------------*/
	public static final String linksFile = "IntLinks.txt";
	public static final String citationsFile = "IntCitations.txt";
	public static int numDocs = 25053;
	

	private int[][] links;
	private int[][] citations;
	
	public IR_PrjPrt2_Auth_Hub_by_Priyank()
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
	
	
	/*-------------------------------------Start of Main Function----------------------------------------------*/
	
	public static void main(String[] args) throws Exception {
		
		//Initializing the numDocs value and creating a new variable l to access the links and citations retrieved above.
		IR_PrjPrt2_Auth_Hub_by_Priyank.numDocs = 25054;
		IR_PrjPrt2_Auth_Hub_by_Priyank l = new IR_PrjPrt2_Auth_Hub_by_Priyank();
		
		//Root Set variable to indicate the size of the root set.
		int rootSet = 10;
		//Array to store the root set retrieved from the TF/IDF similarity step.
		int[] docSet = new int[rootSet];
		
		
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
		
		System.out.println("Generating the system......Please wait.");
		System.out.println("Some information: The size of the root set used in the code is: "+rootSet);
		
		
		System.out.print("query> ");
		
		while(!(str = sc.nextLine()).equals("")){
			System.out.println("Processing Query............");
			
			
		
			
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
				
				if(SoP_TF_IDF != 0)
					MagQueryVec.put(count,SoP_TF_IDF/(disQ * Math.round(Math.sqrt(MagDocVec_TF_IDF.get(count))*100000)/100000));
			}
			
			// calling a sorting function
			HashMap<Integer,Double> FinalSortedList = HashSort(MagQueryVec);
			
			/*-------3. Iterating through the Sorted List and Displaying the Ranked documents------*/
			Set RankedList = FinalSortedList.entrySet();
			loop = RankedList.iterator();
			
			int i = 0;
	/*---------------------------------------------------------End of the code from Project Part 1-----------------------------------------*/		
			
			
			/*
			 *	Loop to store the top K documents set into the docSet, i.e. the root set. 
			*/
			while(loop.hasNext()) {
				Map.Entry finalMap = (Map.Entry)loop.next();
				
				docSet[i] = (int) finalMap.getKey();
				
				if(i == rootSet-1)
					break;
				i++;
			}		
			
			
			//HashMap to store the mapping between the Document ID and the row it belongs to in the Authority and Hub rank array/matrix.
			HashMap<Integer,Integer> AuthorityMapSet = new HashMap<Integer, Integer>();
			
			//HashMap to store the Authority Ranking w.r.t to the document ID
			HashMap<Integer,Double> AuthRankSet = new HashMap<Integer, Double>();
			
			//HashMap to store the Hub Ranking w.r.t to the document ID
			HashMap<Integer,Double> HubRankSet = new HashMap<Integer, Double>();
			
			//variable to store the size of the adjacency matrix.
			int adjMatrixSize = 0;

			//Integer Array to store the base-set documents
			ArrayList<Integer> baseSetDocument = new ArrayList<Integer>();
		    

		/*----------------------------------Loop to generate the base-set by iterating through all the Links and Citations list---------------------------------------*/
			//Here, i generate to unique base-set by checking in the if conditions for its presence.
			for(int Item: docSet){
				if(!baseSetDocument.contains(Item))
					baseSetDocument.add(Item);
				int[] links = l.getLinks(Item);
				for(int link: links){
					if(baseSetDocument.contains(link))
						continue;
					else{
						baseSetDocument.add(link);
					}
				}
				int[] cits = l.getCitations(Item);
				for(int citation: cits){
					if(baseSetDocument.contains(citation))
						continue;
					else{
						baseSetDocument.add(citation);
					}
				}
			}
			
		/*----------------------------------End of base-set generation loop----------------------------------------------------------------*/
			
			//Mapping the base-set documents with its row number in the array. This would be used just to identify the documents after sorting.
			for(int newBaseRow = 0; newBaseRow < baseSetDocument.size(); newBaseRow++){
				AuthorityMapSet.put(newBaseRow, baseSetDocument.get(newBaseRow));
			}

			//initialize the adjMatrixSize variable
			adjMatrixSize = baseSetDocument.size();


			//Creating the Adjacency Matrix by calling the CreateAdjMatrix() function. Definition of this function is provided at the bottom of the program. 
			double[][] adjMatrix = CreateAdjMatrix(baseSetDocument,l);
			
			//Declaring the Transpose of the adjacency Matrix.
			double[][] adjMatrixTrans = new double[adjMatrixSize][adjMatrixSize];
			
			//Initializing all the required arrays before the calculations.
			double[] Auth_curr = CreateInitialMatrix(adjMatrixSize);
			double[] Hub_curr = CreateInitialMatrix(adjMatrixSize);
			double[] Auth_prev = CreateInitialMatrix(adjMatrixSize);
			double[] Hub_prev = CreateInitialMatrix(adjMatrixSize);

			//Transposing the matrix using the user-defined function.
			adjMatrixTrans = transposeMatrix(adjMatrix);

			//initializing the LargestError and the Threshold value.
			double LargestError = 1;
			double Threshold = 0.000000000000001;

		/*---------------------------------------------Convergence loop to calculate the final Authority and the Hubs after going through the power iterations----------------------------*/ 
			while(LargestError > Threshold){

				Auth_prev = equalToArray(Auth_curr);
				Hub_prev = equalToArray(Hub_curr);

				Auth_curr = MatrixMultiply(adjMatrixTrans,Hub_prev);
				Hub_curr = MatrixMultiply(adjMatrix,Auth_curr);

				Auth_curr = NormalizeMatrix(Auth_curr);
				Hub_curr = NormalizeMatrix(Hub_curr);

				LargestError = Math.max(findingLargest(Auth_curr,Auth_prev), findingLargest(Hub_curr,Hub_prev));

			}
		/*----------------------------------------End of Convergence loop------------------------------------------------------------*/
			
			//Storing the values of the Authorities into a HashMap with its resp. document ID. This will be sorted for the final ordering.
			for(int row = 0; row < Auth_curr.length; row++){
				AuthRankSet.put(AuthorityMapSet.get(row), Auth_curr[row]);
			}

			//Final array sorted HashMap
			HashMap<Integer,Double> AuthoritySortedList = HashSort(AuthRankSet);


			/*-------3. Iterating through the Sorted List and Displaying the Ranked documents------*/
			Set AuthRankedList = AuthoritySortedList.entrySet();
			loop = AuthRankedList.iterator();

			i = 0;
			
			System.out.println("Authority Top 10:");
			while(loop.hasNext()) {
				Map.Entry finalMap = (Map.Entry)loop.next();

				//docSet[i] = (int) finalMap.getKey();
				//String d_url = r.document((int) finalMap.getKey()).getFieldable("path").stringValue().replace("%%", "/");
				System.out.println("Doc: "+finalMap.getKey());

				if(i == 9)
					break;
				i++;
			}	
			/*-------------------------------------End of Loop------------------------------------*/
			
			
			//Storing the values of the Hub into a HashMap with its resp. document ID. This will be sorted for the final ordering.
			for(int row = 0; row < Auth_curr.length; row++){
				HubRankSet.put(AuthorityMapSet.get(row), Hub_curr[row]);
			}
			
			//Final Sorted Array
			HashMap<Integer,Double> HubsSortedList = HashSort(HubRankSet);


			/*-------3. Iterating through the Sorted List and Displaying the Ranked documents------*/
			Set HubsRankedList = HubsSortedList.entrySet();
			loop = HubsRankedList.iterator();

			i = 0;

			System.out.println("");
			
			System.out.println("Hubs Top 10:");
			while(loop.hasNext()) {
				Map.Entry finalMap = (Map.Entry)loop.next();

				System.out.println("Doc: "+finalMap.getKey());

				if(i == 9)
					break;
				i++;
			}	
			/*-------------------End of the Loop--------------------------------*/
			
			
			System.out.print("");
			System.out.print("query> ");
			
		}
		
		sc.close();
	}
	
	/*
	 * Creating the Adjacency Matrix
	 *  
	 */
	
	public static double[][] CreateAdjMatrix(ArrayList<Integer> ItemList, IR_PrjPrt2_Auth_Hub_by_Priyank LA){
		double[][] adjMatrixNew = new double[ItemList.size()][ItemList.size()];
		int row = 0, col = 0;
		
		for(row = 0;row < ItemList.size(); row++){
			
			int docNumber = ItemList.get(row);
			int[] links = LA.getLinks(docNumber);
			
			List<Integer> intList = new ArrayList<Integer>(links.length);
			for (int i=0; i<links.length; i++)
			    intList.add(links[i]);
			
			for(col = 0; col < ItemList.size(); col++){
				if(intList.contains(ItemList.get(col)))
					adjMatrixNew[row][col] = 1.0;
				else
					adjMatrixNew[row][col] = 0.0;
			}
		}
		
		
		return adjMatrixNew;
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
	 * Function to transpose the matrix
	 */
	public static double[][] transposeMatrix(double[][] originalMatrix){
		double[][] transposedMatrix = new double[originalMatrix[0].length][originalMatrix.length];
		
		for(int row = 0; row < originalMatrix[0].length; row++)
			for(int col = 0; col < originalMatrix.length; col++)
				transposedMatrix[row][col] = originalMatrix[col][row];
		return transposedMatrix;
	}
	
	/*
	 * Function to perform Matrix Multiplication
	 */
	public static double[] MatrixMultiply(double[][] first, double[] second){
		double[] result = new double[first.length];
		for(int row = 0; row < first.length; row++)
			for(int col = 0; col < first[0].length; col++)
				result[row] = result[row] + (first[row][col]*second[col]);
		
		return result;
	}
	
	/*
	 * Function to Normalize the Authorities and the Hub
	 */
	public static double[] NormalizeMatrix(double[] targetVector){
		double L2Norm = 0;
		for(int row = 0; row < targetVector.length; row++)
			L2Norm = L2Norm + Math.pow(targetVector[row], 2);
		for(int row = 0; row < targetVector.length; row++)
			targetVector[row] = targetVector[row]/(Math.sqrt(L2Norm));
		
		return targetVector;
	}
	
	/*
	 * Function to find the largest value
	 */
	public static double findingLargest(double[] firstVector, double[] secondVector){
		double[] result = new double[firstVector.length];
		double largestValue = 0;
		for(int row = 0; row < firstVector.length; row++){
			result[row] = Math.abs(firstVector[row] - secondVector[row]);
			if(result[row] >= largestValue)
				largestValue = result[row];
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
