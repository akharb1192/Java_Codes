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


class Solution {
    public List<String> topKFrequent(String[] words, int k) {
        HashMap<String,Integer> hmap=new HashMap<>();
        for(String s:words)
        {
            hmap.merge(s,1,Integer::sum);
        }
        PriorityQueue<Map.Entry<String,Integer>> pq=new PriorityQueue<>(k,(a,b)->
            a.getValue()==b.getValue()?a.getKey().compareTo(b.getKey()):b.getValue()-a.getValue()
        );
        for(Map.Entry<String,Integer> m:hmap.entrySet())
        {
            pq.add(m);
        }
        List<String> lis=new ArrayList();
        while(lis.size()<k)
        {
            lis.add(pq.poll().getKey());
        }
        return lis;
    }

}


// Java program to check if there is exist a path between two vertices 
// of a graph. 
import java.io.*; 
import java.util.*; 
import java.util.LinkedList; 

// This class represents a directed graph using adjacency list 
// representation 
class Graph 
{ 
	private int V; // No. of vertices 
	private LinkedList<Integer> adj[]; //Adjacency List 

	//Constructor 
	Graph(int v) 
	{ 
		V = v; 
		adj = new LinkedList[v]; 
		for (int i=0; i<v; ++i) 
			adj[i] = new LinkedList(); 
	} 

	//Function to add an edge into the graph 
	void addEdge(int v,int w) { adj[v].add(w); } 

	//prints BFS traversal from a given source s 
	Boolean isReachable(int s, int d) 
	{ 
		LinkedList<Integer>temp; 

		// Mark all the vertices as not visited(By default set 
		// as false) 
		boolean visited[] = new boolean[V]; 

		// Create a queue for BFS 
		LinkedList<Integer> queue = new LinkedList<Integer>(); 

		// Mark the current node as visited and enqueue it 
		visited[s]=true; 
		queue.add(s); 

		// 'i' will be used to get all adjacent vertices of a vertex 
		Iterator<Integer> i; 
		while (queue.size()!=0) 
		{ 
			// Dequeue a vertex from queue and print it 
			s = queue.poll(); 

			int n; 
			i = adj[s].listIterator(); 

			// Get all adjacent vertices of the dequeued vertex s 
			// If a adjacent has not been visited, then mark it 
			// visited and enqueue it 
			while (i.hasNext()) 
			{ 
				n = i.next(); 

				// If this adjacent node is the destination node, 
				// then return true 
				if (n==d) 
					return true; 

				// Else, continue to do BFS 
				if (!visited[n]) 
				{ 
					visited[n] = true; 
					queue.add(n); 
				} 
			} 
		} 

		// If BFS is complete without visited d 
		return false; 
	} 

	// Driver method 
	public static void main(String args[]) 
	{ 
		// Create a graph given in the above diagram 
		Graph g = new Graph(4); 
		g.addEdge(0, 1); 
		g.addEdge(0, 2); 
		g.addEdge(1, 2); 
		g.addEdge(2, 0); 
		g.addEdge(2, 3); 
		g.addEdge(3, 3); 

		int u = 1; 
		int v = 3; 
		if (g.isReachable(u, v)) 
			System.out.println("There is a path from " + u +" to " + v); 
		else
			System.out.println("There is no path from " + u +" to " + v);; 

		u = 3; 
		v = 1; 
		if (g.isReachable(u, v)) 
			System.out.println("There is a path from " + u +" to " + v); 
		else
			System.out.println("There is no path from " + u +" to " + v);; 
	} 
}


class Solution {
    public String largestNumber(int[] nums) {
        if (nums == null || nums.length == 0) {
            return "";
        }

        // Convert int array to String array, so we can sort later on
        String[] strs = new String[nums.length];
        for (int i = 0; i < nums.length; i++) {
            strs[i] = String.valueOf(nums[i]);
        }

        // Comparator to decide which string should come first in concatenation
        Arrays.sort(strs, (str1, str2) -> (str2 + str1).compareTo(str1 + str2));

        // An extreme edge case that you have only a bunch of 0 in your int array
        if (strs[0].charAt(0) == '0') {
            return "0";
        }
        StringBuilder sb = new StringBuilder();
        for (String s : strs) {
            sb.append(s);
        }
        return sb.toString();
    }
}

class Solution {
    public boolean wordBreak(String s, List<String> wordDict) {
        Set<String> wordDictSet = new HashSet(wordDict);
        boolean[] dp = new boolean[s.length() + 1];
        dp[0] = true;
        for (int i = 1; i <= s.length(); i++) {
            for (int j = 0; j < i; j++) {
                if (dp[j] && wordDictSet.contains(s.substring(j, i))) {
                    dp[i] = true;
                    break;
                }
            }
        }
        return dp[s.length()];
    }
}



class Solution {
    int maxDownVal;
    public int maxPathSum(TreeNode root) {
        maxDownVal = Integer.MIN_VALUE;
        updateVal(root);
        return maxDownVal;
    }
    
    private int updateVal(TreeNode root) {
        if (root == null) {
            return 0;
        }
        
        int left = Math.max(0, updateVal(root.left));
        int right = Math.max(0, updateVal(root.right));
        
        maxDownVal = Math.max(maxDownVal, left+right+root.val);
        
        return Math.max(left, right) + root.val;
    }
}



// Java program to find next greater 
// number with same set of digits. 
import java.util.Arrays; 

public class nextGreater 
{ 
	// Utility function to swap two digit 
	static void swap(char ar[], int i, int j) 
	{ 
		char temp = ar[i]; 
		ar[i] = ar[j]; 
		ar[j] = temp; 
	} 

	// Given a number as a char array number[], 
	// this function finds the next greater number. 
	// It modifies the same array to store the result 
	static void findNext(char ar[], int n) 
	{ 
		int i; 
		
		// I) Start from the right most digit 
		// and find the first digit that is smaller 
		// than the digit next to it. 
		for (i = n - 1; i > 0; i--) 
		{ 
			if (ar[i] > ar[i - 1]) { 
				break; 
			} 
		} 
		
		// If no such digit is found, then all 
		// digits are in descending order means 
		// there cannot be a greater number with 
		// same set of digits 
		if (i == 0) 
		{ 
			System.out.println("Not possible"); 
		} 
		else
		{ 
			int x = ar[i - 1], min = i; 
			
			// II) Find the smallest digit on right 
			// side of (i-1)'th digit that is greater 
			// than number[i-1] 
			for (int j = i + 1; j < n; j++) 
			{ 
				if (ar[j] > x && ar[j] < ar[min]) 
				{ 
					min = j; 
				} 
			} 

			// III) Swap the above found smallest 
			// digit with number[i-1] 
			swap(ar, i - 1, min); 

			// IV) Sort the digits after (i-1) 
			// in ascending order 
			Arrays.sort(ar, i, n); 
			System.out.print("Next number with same" + 
									" set of digits is "); 
			for (i = 0; i < n; i++) 
				System.out.print(ar[i]); 
		} 
	} 

	public static void main(String[] args) 
	{ 
		char digits[] = { '5','3','4','9','7','6' }; 
		int n = digits.length; 
		findNext(digits, n); 
	} 
}



class Solution {
    public ListNode addTwoNumbers(ListNode m, ListNode n) {
        int carry = 0;
        ListNode dummy = new ListNode(0);
        ListNode tail = dummy;
        while (m != null || n!= null || carry != 0) {
            int value = carry;
            if (m != null) {
                value += m.val;
                m = m.next;
            }
            if (n != null) {
                value += n.val;
                n = n.next;
            }
            int digit = value % 10;
            carry = value / 10;
            tail.next = new ListNode(digit);
            tail = tail.next;
        }
        return dummy.next;
    }
}


class Solution {
    public ListNode addTwoNumbers(ListNode m, ListNode n) {
        if (m == null) {
            return n;
        } else if (n == null) {
            return m;
        }

        Deque<Integer> deque1 = new ArrayDeque(); 
        Deque<Integer> deque2 = new ArrayDeque(); 

        while (m != null) {
            deque1.push(m.val);
            m = m.next;
        }
        while (n != null) {
            deque2.push(n.val);
            n = n.next;
        }

        int carry = 0;
        ListNode dummy = new ListNode(0);   
        while (!deque1.isEmpty() || !deque2.isEmpty() || carry != 0) {
            int value = carry;
            if (!deque1.isEmpty()) {
                value += deque1.pop();
            }
            if (!deque2.isEmpty()) {
                value += deque2.pop();
            }
            int digit = value % 10;
            carry = value / 10;
            ListNode head = new ListNode(digit);
            head.next = dummy.next;
            dummy.next = head;
        }
        return dummy.next;
    }
}


class Solution {
    public boolean checkInclusion(String s1, String s2) {
        // Init the table that keeps track of char in s1
        int[] table = new int[26];
        int count = s1.length();
        for(int i = 0; i < s1.length(); i++) {
            table[s1.charAt(i) - 'a']++;
        }
        
        // Fixed length window is [begin, begin + s1.length())
        int begin = 0, end = s1.length();
        for(int i = 0; i < end && i < s2.length(); i++) {
            if(table[s2.charAt(i) - 'a']-- > 0) {
                count--;
            }
        }
        if(count == 0) return true;
        
        // Move the window to the right
        while(end < s2.length()) {
            if(table[s2.charAt(end++) - 'a']-- > 0) {
                count--;
            }
            
            if(table[s2.charAt(begin++) - 'a']++ >= 0) {
                count++;
            }
            
            if(count == 0) return true;
            
        }
        
        return false;
    }
}


public class Solution {
    public String longestPalindrome(String s) {
        char[] cs = s.toCharArray();
        int n = cs.length;
        
        int maxLen = 0;
        int maxI = 0;
        
        for(int i = 0; i < n; i++){
            for(int offset = 0; offset <= 1; offset++){
                int left = i;
                int right = i + offset;
                
                while(left >= 0 && right < n && cs[left] == cs[right]){
                    left--;
                    right++;
                }
                left++;
                right--;
                int len = right - left + 1;
                if(len > maxLen){
                    maxLen = len;
                    maxI = left;
                }
            }
        }
        return s.substring(maxI, maxI + maxLen);
    }
}


class Solution {
    public String minWindow(String s, String t) {
        Map<Character, Integer> map = new HashMap<>();
        for (char c : t.toCharArray()) {
            map.put(c, (map.getOrDefault(c, 0) + 1));
        }

        int start = 0;
        int end = 0;
        int count = map.size();
        int minLen = Integer.MAX_VALUE;
        String ans = "";

        while (end < s.length()) {
            if (map.containsKey(s.charAt(end))) {
                map.put(s.charAt(end), map.get(s.charAt(end)) - 1);
                if (map.get(s.charAt(end)) == 0) {
                    count--;
                }
            }

            end++;

            while (count == 0) {
                if (end - start < minLen) {
                    minLen = end - start;
                    ans = s.substring(start, end);
                }

                if (map.containsKey(s.charAt(start))) {
                    map.put(s.charAt(start), map.get(s.charAt(start)) + 1);
                    if (map.get(s.charAt(start)) > 0) {
                        count++;
                    }
                }

                start++;
            }
        }

        return ans;
    }
}


public TreeNode lowestCommonAncestor(TreeNode root, TreeNode p, TreeNode q) {
    if(root==null)
        return null;
 
    if(root==p || root==q)
        return root;
 
    TreeNode l = lowestCommonAncestor(root.left, p, q);
    TreeNode r = lowestCommonAncestor(root.right, p, q);
 
    if(l!=null&&r!=null){
        return root;
    }else if(l==null&&r==null){
        return null;
    }else{
        return l==null?r:l;
    }
}

// Java program to search 
// a word in a 2D grid 
import java.io.*; 
import java.util.*; 

class GFG { 

	// Rows and columns in the given grid 
	static int R, C; 

	// For searching in all 8 direction 
	static int[] x = { -1, -1, -1, 0, 0, 1, 1, 1 }; 
	static int[] y = { -1, 0, 1, -1, 1, -1, 0, 1 }; 

	// This function searches in all 
	// 8-direction from point 
	// (row, col) in grid[][] 
	static boolean search2D(char[][] grid, int row, 
							int col, String word) 
	{ 
		// If first character of word 
		// doesn't match with 
		// given starting point in grid. 
		if (grid[row][col] != word.charAt(0)) 
			return false; 

		int len = word.length(); 

		// Search word in all 8 directions 
		// starting from (row, col) 
		for (int dir = 0; dir < 8; dir++) { 
			// Initialize starting point 
			// for current direction 
			int k, rd = row + x[dir], cd = col + y[dir]; 

			// First character is already checked, 
			// match remaining characters 
			for (k = 1; k < len; k++) { 
				// If out of bound break 
				if (rd >= R || rd < 0 || cd >= C || cd < 0) 
					break; 

				// If not matched, break 
				if (grid[rd][cd] != word.charAt(k)) 
					break; 

				// Moving in particular direction 
				rd += x[dir]; 
				cd += y[dir]; 
			} 

			// If all character matched, 
			// then value of must 
			// be equal to length of word 
			if (k == len) 
				return true; 
		} 
		return false; 
	} 

	// Searches given word in a given 
	// matrix in all 8 directions 
	static void patternSearch( 
		char[][] grid, 
		String word) 
	{ 
		// Consider every point as starting 
		// point and search given word 
		for (int row = 0; row < R; row++) { 
			for (int col = 0; col < C; col++) { 
				if (search2D(grid, row, col, word)) 
					System.out.println( 
						"pattern found at " + row + ", " + col); 
			} 
		} 
	} 

	// Driver code 
	public static void main(String args[]) 
	{ 
		R = 3; 
		C = 13; 
		char[][] grid = { { 'G', 'E', 'E', 'K', 'S', 'F', 'O', 'R', 'G', 'E', 'E', 'K', 'S' }, 
						{ 'G', 'E', 'E', 'K', 'S', 'Q', 'U', 'I', 'Z', 'G', 'E', 'E', 'K' }, 
						{ 'I', 'D', 'E', 'Q', 'A', 'P', 'R', 'A', 'C', 'T', 'I', 'C', 'E' } }; 
		patternSearch(grid, "GEEKS"); 
		System.out.println(); 
		patternSearch(grid, "EEE"); 
	} 
} 



/* Iterative Java program to find sum of data of all leaves 
of a binary tree on same level and then multiply sums 
obtained of all levels. */

/* importing the necessary class */
import java.util.LinkedList; 
import java.util.Queue; 
import java.util.Stack; 

/* Class containing left and right child of current 
node and key value*/
class Node { 

	int data; 
	Node left, right; 

	public Node(int item) 
	{ 
		data = item; 
		left = right = null; 
	} 
} 

class BinaryTree { 

	Node root; 

	// helper function to check if a Node is leaf of tree 
	boolean isLeaf(Node node) 
	{ 
		return ((node.left == null) && (node.right == null)); 
	} 

	/* Calculate sum of all leaf Nodes at each level and returns 
	multiplication of sums */
	int sumAndMultiplyLevelData() 
	{ 
		return sumAndMultiplyLevelData(root); 
	} 
	int sumAndMultiplyLevelData(Node node) 
	{ 
		// Tree is empty 
		if (node == null) { 
			return 0; 
		} 

		int mul = 1; /* To store result */

		// Create an empty queue for level order tarversal 
		LinkedList<Node> q = new LinkedList<Node>(); 

		// Enqueue Root and initialize height 
		q.add(node); 

		// Do level order traversal of tree 
		while (true) { 

			// NodeCount (queue size) indicates number of Nodes 
			// at current lelvel. 
			int NodeCount = q.size(); 

			// If there are no Nodes at current level, we are done 
			if (NodeCount == 0) { 
				break; 
			} 

			// Initialize leaf sum for current level 
			int levelSum = 0; 

			// A boolean variable to indicate if found a leaf 
			// Node at current level or not 
			boolean leafFound = false; 

			// Dequeue all Nodes of current level and Enqueue all 
			// Nodes of next level 
			while (NodeCount > 0) { 
				Node node1; 
				node1 = q.poll(); 

				/* if Node is a leaf, update sum at the level */
				if (isLeaf(node1)) { 
					leafFound = true; 
					levelSum += node1.data; 
				} 

				// Add children of Node 
				if (node1.left != null) { 
					q.add(node1.left); 
				} 
				if (node1.right != null) { 
					q.add(node1.right); 
				} 
				NodeCount--; 
			} 

			// If we found at least one leaf, we multiply 
			// result with level sum. 
			if (leafFound) { 
				mul *= levelSum; 
			} 
		} 

		return mul; // Return result 
	} 

	public static void main(String args[]) 
	{ 

		/* creating a binary tree and entering 
		the nodes */
		BinaryTree tree = new BinaryTree(); 
		tree.root = new Node(2); 
		tree.root.left = new Node(7); 
		tree.root.right = new Node(5); 
		tree.root.left.left = new Node(8); 
		tree.root.left.right = new Node(6); 
		tree.root.left.right.left = new Node(1); 
		tree.root.left.right.right = new Node(11); 
		tree.root.right.right = new Node(9); 
		tree.root.right.right.left = new Node(4); 
		tree.root.right.right.right = new Node(10); 
		System.out.println("The final product value : "
						+ tree.sumAndMultiplyLevelData()); 
	} 
}


import java.util.*; 
import java.io.*; 
class Node { 
	int data; 
	Node left, right, nextRight; 

	Node(int item) 
	{ 
		data = item; 
		left = right = nextRight = null; 
	} 
} 

public class BinaryTree { 
	Node root; 
	void connect(Node p) 
	{ 
		// initialize queue to hold nodes at same level 
		Queue<Node> q = new LinkedList<>(); 

		q.add(root); // adding nodes to tehe queue 

		Node temp = null; // initializing prev to null 
		while (!q.isEmpty()) { 
			int n = q.size(); 
			for (int i = 0; i < n; i++) { 
				Node prev = temp; 
				temp = q.poll(); 

				// i > 0 because when i is 0 prev points 
				// the last node of previous level, 
				// so we skip it 
				if (i > 0) 
					prev.nextRight = temp; 

				if (temp.left != null) 
					q.add(temp.left); 

				if (temp.right != null) 
					q.add(temp.right); 
			} 

			// pointing last node of the nth level to null 
			temp.nextRight = null; 
		} 
	} 

	// Driver program to test above functions 
	public static void main(String args[]) 
	{ 
		BinaryTree tree = new BinaryTree(); 

		/* Constructed binary tree is 
			10 
			/ \ 
		8	 2 
		/ 
		3 
		*/
		tree.root = new Node(10); 
		tree.root.left = new Node(8); 
		tree.root.right = new Node(2); 
		tree.root.left.left = new Node(3); 

		// Populates nextRight pointer in all nodes 
		tree.connect(tree.root); 

		// Let us check the values of nextRight pointers 
		System.out.println("Following are populated nextRight pointers in "
						+ "the tree"
						+ "(-1 is printed if there is no nextRight)"); 
		int a = tree.root.nextRight != null ? tree.root.nextRight.data : -1; 
		System.out.println("nextRight of " + tree.root.data + " is "
						+ a); 
		int b = tree.root.left.nextRight != null ? tree.root.left.nextRight.data : -1; 
		System.out.println("nextRight of " + tree.root.left.data + " is "
						+ b); 
		int c = tree.root.right.nextRight != null ? tree.root.right.nextRight.data : -1; 
		System.out.println("nextRight of " + tree.root.right.data + " is "
						+ c); 
		int d = tree.root.left.left.nextRight != null ? tree.root.left.left.nextRight.data : -1; 
		System.out.println("nextRight of " + tree.root.left.left.data + " is "
						+ d); 
	} 
} 



public class Solution {
    public boolean isBalanced(TreeNode root) {
        if(root==null)
        {
            return true;
        }
        return isBalanced(root.left) && isBalanced(root.right) && Math.abs(depth(root.left)-depth(root.right))<=1;
    }
    
    public int depth(TreeNode root) {
       if(root==null)
       {
           return 0;
       }
        return 1+Math.max(depth(root.left),depth(root.right));
    }
}


// Java program to construct a tree using inorder and preorder traversal 

/* A binary tree node has data, pointer to left child 
and a pointer to right child */
class Node { 
	char data; 
	Node left, right; 

	Node(char item) 
	{ 
		data = item; 
		left = right = null; 
	} 
} 

class BinaryTree { 
	Node root; 
	static int preIndex = 0; 

	/* Recursive function to construct binary of size len from 
	Inorder traversal in[] and Preorder traversal pre[]. 
	Initial values of inStrt and inEnd should be 0 and len -1. 
	The function doesn't do any error checking for cases where 
	inorder and preorder do not form a tree */
	Node buildTree(char in[], char pre[], int inStrt, int inEnd) 
	{ 
		if (inStrt > inEnd) 
			return null; 

		/* Pick current node from Preorder traversal using preIndex 
		and increment preIndex */
		Node tNode = new Node(pre[preIndex++]); 

		/* If this node has no children then return */
		if (inStrt == inEnd) 
			return tNode; 

		/* Else find the index of this node in Inorder traversal */
		int inIndex = search(in, inStrt, inEnd, tNode.data); 

		/* Using index in Inorder traversal, construct left and 
		right subtress */
		tNode.left = buildTree(in, pre, inStrt, inIndex - 1); 
		tNode.right = buildTree(in, pre, inIndex + 1, inEnd); 

		return tNode; 
	} 

	/* UTILITY FUNCTIONS */

	/* Function to find index of value in arr[start...end] 
	The function assumes that value is present in in[] */
	int search(char arr[], int strt, int end, char value) 
	{ 
		int i; 
		for (i = strt; i <= end; i++) { 
			if (arr[i] == value) 
				return i; 
		} 
		return i; 
	} 

	/* This funtcion is here just to test buildTree() */
	void printInorder(Node node) 
	{ 
		if (node == null) 
			return; 

		/* first recur on left child */
		printInorder(node.left); 

		/* then print the data of node */
		System.out.print(node.data + " "); 

		/* now recur on right child */
		printInorder(node.right); 
	} 

	// driver program to test above functions 
	public static void main(String args[]) 
	{ 
		BinaryTree tree = new BinaryTree(); 
		char in[] = new char[] { 'D', 'B', 'E', 'A', 'F', 'C' }; 
		char pre[] = new char[] { 'A', 'B', 'D', 'E', 'C', 'F' }; 
		int len = in.length; 
		Node root = tree.buildTree(in, pre, 0, len - 1); 

		// building the tree by printing inorder traversal 
		System.out.println("Inorder traversal of constructed tree is : "); 
		tree.printInorder(root); 
	} 
}



// Java program to print all nodes at a distance k from given node 

// A binary tree node 
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
	/* Recursive function to print all the nodes at distance k in 
	tree (or subtree) rooted with given root. */

	void printkdistanceNodeDown(Node node, int k) 
	{ 
		// Base Case 
		if (node == null || k < 0) 
			return; 

		// If we reach a k distant node, print it 
		if (k == 0) 
		{ 
			System.out.print(node.data); 
			System.out.println(""); 
			return; 
		} 

		// Recur for left and right subtrees 
		printkdistanceNodeDown(node.left, k - 1); 
		printkdistanceNodeDown(node.right, k - 1); 
	} 

	// Prints all nodes at distance k from a given target node. 
	// The k distant nodes may be upward or downward.This function 
	// Returns distance of root from target node, it returns -1 
	// if target node is not present in tree rooted with root. 
	int printkdistanceNode(Node node, Node target, int k) 
	{ 
		// Base Case 1: If tree is empty, return -1 
		if (node == null) 
			return -1; 

		// If target is same as root. Use the downward function 
		// to print all nodes at distance k in subtree rooted with 
		// target or root 
		if (node == target) 
		{ 
			printkdistanceNodeDown(node, k); 
			return 0; 
		} 

		// Recur for left subtree 
		int dl = printkdistanceNode(node.left, target, k); 

		// Check if target node was found in left subtree 
		if (dl != -1) 
		{ 
			
			// If root is at distance k from target, print root 
			// Note that dl is Distance of root's left child from 
			// target 
			if (dl + 1 == k) 
			{ 
				System.out.print(node.data); 
				System.out.println(""); 
			} 
			
			// Else go to right subtree and print all k-dl-2 distant nodes 
			// Note that the right child is 2 edges away from left child 
			else
				printkdistanceNodeDown(node.right, k - dl - 2); 

			// Add 1 to the distance and return value for parent calls 
			return 1 + dl; 
		} 

		// MIRROR OF ABOVE CODE FOR RIGHT SUBTREE 
		// Note that we reach here only when node was not found in left 
		// subtree 
		int dr = printkdistanceNode(node.right, target, k); 
		if (dr != -1) 
		{ 
			if (dr + 1 == k) 
			{ 
				System.out.print(node.data); 
				System.out.println(""); 
			} 
			else
				printkdistanceNodeDown(node.left, k - dr - 2); 
			return 1 + dr; 
		} 

		// If target was neither present in left nor in right subtree 
		return -1; 
	} 

	// Driver program to test the above functions 
	public static void main(String args[]) 
	{ 
		BinaryTree tree = new BinaryTree(); 

		/* Let us construct the tree shown in above diagram */
		tree.root = new Node(20); 
		tree.root.left = new Node(8); 
		tree.root.right = new Node(22); 
		tree.root.left.left = new Node(4); 
		tree.root.left.right = new Node(12); 
		tree.root.left.right.left = new Node(10); 
		tree.root.left.right.right = new Node(14); 
		Node target = tree.root.left.right; 
		tree.printkdistanceNode(tree.root, target, 2); 
	} 
} 



class Solution {
    public List<List<Integer>> subsets(int[] array) {
        if (array == null || array.length == 0) {
            return new ArrayList();
        }
        List<List<Integer>> solutions = new ArrayList();
        makeSubsets(array, 0, solutions, new ArrayList());
        return solutions;
    }

    private void makeSubsets(int[] array, int i, List<List<Integer>> solutions, List<Integer> list)     {
        if (i == array.length) {
            solutions.add(new ArrayList(list));
            return;
        }

        // don't use array[i]
        makeSubsets(array, i + 1, solutions, list);

        // use array[i]
        list.add(array[i]);
        makeSubsets(array, i + 1, solutions, list);
        list.remove(list.size() - 1);
    }
}



class Solution {
    public List<List<Integer>> permute(int[] array) {
        if (array == null || array.length == 0) {
            return new ArrayList();
        }
        List<List<Integer>> solutions = new ArrayList();
        permute(array, 0, new boolean[array.length], solutions, new ArrayList());
        return solutions;
    }

    private void permute(int[] array, int index, boolean[] used, List<List<Integer>> solutions, List<Integer> list) {
        if (index == array.length) {
            solutions.add(new ArrayList(list));
            return;
        }
        for (int i = 0; i < array.length; i++) {
            if (used[i] == false) {
                list.add(array[i]);
                used[i] = true;
                permute(array, index + 1, used, solutions, list);
                used[i] = false;
                list.remove(list.size() - 1);
            }
        }
    }
}



static void rotateMatrix( 
        int N, int mat[][]) 
    { 
        // Consider all squares one by one 
        for (int x = 0; x < N / 2; x++) { 
            // Consider elements in group 
            // of 4 in current square 
            for (int y = x; y < N - x - 1; y++) { 
                // Store current cell in 
                // temp variable 
                int temp = mat[x][y]; 
  
                // Move values from right to top 
                mat[x][y] = mat[y][N - 1 - x]; 
  
                // Move values from bottom to right 
                mat[y][N - 1 - x] 
                    = mat[N - 1 - x][N - 1 - y]; 
  
                // Move values from left to bottom 
                mat[N - 1 - x][N - 1 - y] = mat[N - 1 - y][x]; 
  
                // Assign temp to left 
                mat[N - 1 - y][x] = temp; 
            } 
        } 
    }





public class Solution {
    public int maxCoins(int[] nums) {
        // Extend list with head and tail (both are 1), index starts at 1
        int array[] = new int[nums.length + 2];
        array[0] = 1;
        array[array.length-1] = 1;
        for (int i = 0; i < nums.length; i++) {
            array[i+1] = nums[i];
        }

        // Initialize DP arrays, 1 index based
        int dp[][] = new int[array.length][array.length]; //dp[i][j] stands for max coins at i step, burst j
        for (int i =0; i < array.length; i++) {
            for (int j = 0; j < array.length; j++) {
                dp[i][j] = 0;
            }
        }

        for (int i=1; i< array.length-1; i++) {
            for (int j=i; j >=1; j--) {
                // k as last
                for (int k=j; k <= i; k++) {
                    dp[j][i] = Math.max(array[j-1]*array[k]*array[i+1] + dp[j][k-1] + dp[k+1][i], dp[j][i]);
                }
            }
        }

        return dp[1][array.length-2];
    }
}



import java.util.*;
import java.util.stream.*;
import java.util.Map.Entry;

public class Main
{
	public static void main(String[] args) {
		String str="aaannkkiitttt";
		Map<Character,Integer> hmap=new HashMap<Character,Integer>();
		for(char a:str.toCharArray())
		{
		    hmap.merge(a,1,Integer::sum);
		}
		LinkedHashMap<Character, Integer> sortedMap = 
        hmap.entrySet().stream().
        sorted(Entry.comparingByValue(Collections.reverseOrder())).
        collect(Collectors.toMap(Entry::getKey, Entry::getValue,
                             (e1, e2) -> e1, LinkedHashMap::new));
		System.out.println(sortedMap);
		
	}
}


//preOrder to postOrder and Inorder
