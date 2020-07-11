public class Solution {
    public int minDistance(String word1, String word2) {
        if (word1 == null || word2 == null){
            return -1;
        }
        int[][] distance = new int[word1.length() + 1][word2.length() + 1];
        
        for (int i = 0; i < word2.length() + 1; i++){
            distance[0][i] = i;
        }
        for (int j = 0; j < word1.length() + 1; j++){
            distance[j][0] = j;
        }
        for (int i = 1; i < word1.length() + 1; i++){
            for (int j = 1; j < word2.length() + 1; j++){
                if (word1.charAt(i - 1) == word2.charAt(j - 1)){
                    distance[i][j] = distance[i - 1][j - 1];
                }
                else{
                    distance[i][j] = Math.min(Math.min(distance[i - 1][j - 1], distance[i - 1][j]), distance[i][j - 1]) + 1;
                }
            }
        }
        
        return distance[word1.length()][word2.length()];
    }
}

class Solution {
public int lengthOfLIS(int[] nums) {
    if(nums==null || nums.length==0)
        return 0;
    int max[]=new int[nums.length];
    Arrays.fill(max,1);
    int result=1;
    for(int i=0;i<nums.length;i++)
    {
        for(int j=0;j<i;j++)
        {
            if(nums[i]>nums[j])
            {
                max[i]=Math.max(max[i],max[j]+1);
            }
        }
        result=Math.max(max[i],result);
    }
    return result;
}}


public class Solution {
	public int longestCommonSubsequence(String text1, String text2) {
		int length1 = text1.length();
		int length2 = text2.length();

		int[][] maxLengths = new int[length1 + 1][length2 + 1];
		for (int i = 0; i <= length1; i++) {
			for (int j = 0; j <= length2; j++) {
				if (i != 0 && j != 0) {
					if (text1.charAt(i - 1) == text2.charAt(j - 1)) {
						maxLengths[i][j] = maxLengths[i - 1][j - 1] + 1;
					} else {
						maxLengths[i][j] = Math.max(maxLengths[i - 1][j], maxLengths[i][j - 1]);
					}
				}
			}
		}

		return maxLengths[length1][length2];
	}
}

// Java program to find nth ugly number 
class GFG { 
	
	/*This function divides a by greatest 
	divisible power of b*/
	static int maxDivide(int a, int b) 
	{ 
		while(a % b == 0) 
			a = a/b; 
		return a; 
	} 
	
	/* Function to check if a number 
	is ugly or not */
	static int isUgly(int no) 
	{ 
		no = maxDivide(no, 2); 
		no = maxDivide(no, 3); 
		no = maxDivide(no, 5); 
		
		return (no == 1)? 1 : 0; 
	} 
	
	/* Function to get the nth ugly 
	number*/
	static int getNthUglyNo(int n) 
	{ 
		int i = 1; 
		
		// ugly number count 
		int count = 1; 
		
		// check for all integers 
		// until count becomes n 
		while(n > count) 
		{ 
			i++; 
			if(isUgly(i) == 1) 
				count++; 
		} 
		return i; 
	} 
	
	/* Driver program to test above 
	functions */
	public static void main(String args[]) 
	{ 
		int no = getNthUglyNo(150); 
		System.out.println("150th ugly "
					+ "no. is "+ no); 
	} 
}

// Java program to print all paths with sum k. 
import java.util.*; 

class GFG 
{ 
	
//utility function to print contents of 
//a vector from index i to it's end 
static void printVector( Vector<Integer> v, int i) 
{ 
	for (int j = i; j < v.size(); j++) 
		System.out.print( v.get(j) + " "); 
		System.out.println(); 
} 

// binary tree node 
static class Node 
{ 
	int data; 
	Node left,right; 
	Node(int x) 
	{ 
		data = x; 
		left = right = null; 
	} 
}; 
static Vector<Integer> path = new Vector<Integer>(); 

// This function prints all paths that have sum k 
static void printKPathUtil(Node root, int k) 
{ 
	// empty node 
	if (root == null) 
		return; 

	// add current node to the path 
	path.add(root.data); 

	// check if there's any k sum path 
	// in the left sub-tree. 
	printKPathUtil(root.left, k); 

	// check if there's any k sum path 
	// in the right sub-tree. 
	printKPathUtil(root.right, k); 

	// check if there's any k sum path that 
	// terminates at this node 
	// Traverse the entire path as 
	// there can be negative elements too 
	int f = 0; 
	for (int j = path.size() - 1; j >= 0; j--) 
	{ 
		f += path.get(j); 

		// If path sum is k, print the path 
		if (f == k) 
			printVector(path, j); 
	} 

	// Remove the current element from the path 
	path.remove(path.size() - 1); 
} 

// A wrapper over printKPathUtil() 
static void printKPath(Node root, int k) 
{ 
	path = new Vector<Integer>(); 
	printKPathUtil(root, k); 
} 

// Driver code 
public static void main(String args[]) 
{ 
	Node root = new Node(1); 
	root.left = new Node(3); 
	root.left.left = new Node(2); 
	root.left.right = new Node(1); 
	root.left.right.left = new Node(1); 
	root.right = new Node(-1); 
	root.right.left = new Node(4); 
	root.right.left.left = new Node(1); 
	root.right.left.right = new Node(2); 
	root.right.right = new Node(5); 
	root.right.right.right = new Node(2); 

	int k = 5; 
	printKPath(root, k); 
} 
}


// Java program to find nth number that contains 
// the digit k or divisible by k. 
import java.io.*; 

class GFG 
{ 
	// Function for checking if digit k 
	// is in n or not 
	public static boolean checkdigit(int n, int k) 
	{ 
		while (n != 0) 
		{ 
			// finding remainder 
			int rem = n % 10; 
	
			// if digit found 
			if (rem == k) 
				return true; 
	
			n = n / 10; 
		} 

		return false; 
	} 

	// Function for finding nth number 
	public static int findNthNumber(int n, int k) 
	{ 
		// since k is the first which satisfy th 
	// criteria, so consider it in count making count = 1 
	// and starting from i = k + 1 
		for (int i = k + 1, count = 1; count < n; i++) 
		{ 
		// checking that the number contain 
		// k digit or divisible by k 
		if (checkdigit(i, k) || (i % k == 0)) 
			count++; 

		if (count == n) 
		return i; 
		} 

	return -1; 
	} 

	// Driver code 
	public static void main (String[] args) 
	{ 
		int n = 10, k = 2; 
		System.out.println(findNthNumber(n, k)); 
		
	} 
} 


class Solution {
    public int[] topKFrequent(int[] nums, int k) {
        if (nums == null || nums.length == 0) {
            return new int[0];
        }

        Map<Integer, Integer> map = new HashMap();
        for (int num : nums) {
            map.merge(num, 1, Integer::sum);
        }
        List<List<Integer>> buckets = new ArrayList(nums.length + 1); 
        for (int i = 0; i < nums.length + 1; i++) {
            buckets.add(new ArrayList<Integer>());
        }

        for (Integer num : map.keySet()) {
            int frequency = map.get(num);
            List<Integer> bucket = buckets.get(frequency);
            bucket.add(num);
        }

        List<Integer> solution = new ArrayList();
        for (int i = buckets.size() - 1; i >= 0 && solution.size() < k; i--) {
            List<Integer> bucket = buckets.get(i);
            solution.addAll(bucket);
        }
        List<Integer> ans= solution.subList(0, k); 
        int[] anss=new int[ans.size()];
        for(int i=0;i<ans.size();i++)
        {
            anss[i]=ans.get(i);
        }
        return anss;
    }
}


class Solution {
    public List<List<String>> groupAnagrams(String[] array) {
        if (array == null || array.length == 0) {
            return new ArrayList();
        }

        Map<String, List<String>> map = new HashMap();

        for (String str : array) {
            String key = sortChars(str);
            map.putIfAbsent(key, new ArrayList<String>());
            List<String> anagrams = map.get(key);
            anagrams.add(str);
        }

        return new ArrayList(map.values());
    }

    private String sortChars(String str) {
        char[] content = str.toCharArray(); // Strings are immutable, so we convert to char[]
        Arrays.sort(content);
        return new String(content);
    }
}


class Node 
{ 
	int data; 
	Node left, right; 

	Node(int item) 
	{ 
		data = item; 
		left = right = null; 
	} 
} 

class BinaryTree 
{ 
	Node root; 

	// Convert a given tree to a tree where every node contains sum of 
	// values of nodes in left and right subtrees in the original tree 
	int toSumTree(Node node) 
	{ 
		// Base case 
		if (node == null) 
			return 0; 

		// Store the old value 
		int old_val = node.data; 

		// Recursively call for left and right subtrees and store the sum 
		// as new value of this node 
		node.data = toSumTree(node.left) + toSumTree(node.right); 

		// Return the sum of values of nodes in left and right subtrees 
		// and old_value of this node 
		return node.data + old_val; 
	} 

	// A utility function to print inorder traversal of a Binary Tree 
	void printInorder(Node node) 
	{ 
		if (node == null) 
			return; 
		printInorder(node.left); 
		System.out.print(node.data + " "); 
		printInorder(node.right); 
	} 

	/* Driver function to test above functions */
	public static void main(String args[]) 
	{ 
		BinaryTree tree = new BinaryTree(); 

		/* Constructing tree given in the above figure */
		tree.root = new Node(10); 
		tree.root.left = new Node(-2); 
		tree.root.right = new Node(6); 
		tree.root.left.left = new Node(8); 
		tree.root.left.right = new Node(-4); 
		tree.root.right.left = new Node(7); 
		tree.root.right.right = new Node(5); 

		tree.toSumTree(tree.root); 

		// Print inorder traversal of the converted tree to test result 
		// of toSumTree() 
		System.out.println("Inorder Traversal of the resultant tree is:"); 
		tree.printInorder(tree.root); 
	} 
}



import java.util.ArrayList; 
import java.util.LinkedList; 
import java.util.List; 

class Graph { 
	
	private final int V; 
	private final List<List<Integer>> adj; 

	public Graph(int V) 
	{ 
		this.V = V; 
		adj = new ArrayList<>(V); 
		
		for (int i = 0; i < V; i++) 
			adj.add(new LinkedList<>()); 
	} 
	
	// This function is a variation of DFSUtil() in 
	// https://www.geeksforgeeks.org/archives/18212 
	private boolean isCyclicUtil(int i, boolean[] visited, 
									boolean[] recStack) 
	{ 
		
		// Mark the current node as visited and 
		// part of recursion stack 
		if (recStack[i]) 
			return true; 

		if (visited[i]) 
			return false; 
			
		visited[i] = true; 

		recStack[i] = true; 
		List<Integer> children = adj.get(i); 
		
		for (Integer c: children) 
			if (isCyclicUtil(c, visited, recStack)) 
				return true; 
				
		recStack[i] = false; 

		return false; 
	} 

	private void addEdge(int source, int dest) { 
		adj.get(source).add(dest); 
	} 

	// Returns true if the graph contains a 
	// cycle, else false. 
	// This function is a variation of DFS() in 
	// https://www.geeksforgeeks.org/archives/18212 
	private boolean isCyclic() 
	{ 
		
		// Mark all the vertices as not visited and 
		// not part of recursion stack 
		boolean[] visited = new boolean[V]; 
		boolean[] recStack = new boolean[V]; 
		
		
		// Call the recursive helper function to 
		// detect cycle in different DFS trees 
		for (int i = 0; i < V; i++) 
			if (isCyclicUtil(i, visited, recStack)) 
				return true; 

		return false; 
	} 

	// Driver code 
	public static void main(String[] args) 
	{ 
		Graph graph = new Graph(4); 
		graph.addEdge(0, 1); 
		graph.addEdge(0, 2); 
		graph.addEdge(1, 2); 
		graph.addEdge(2, 0); 
		graph.addEdge(2, 3); 
		graph.addEdge(3, 3); 
		
		if(graph.isCyclic()) 
			System.out.println("Graph contains cycle"); 
		else
			System.out.println("Graph doesn't "
									+ "contain cycle"); 
	} 
}
