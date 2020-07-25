//Permutations of a string

import java.util.*;
public class Main
{
    public static void permutation(String str,int i,int j)
    {
        if(i==j)
        {
            System.out.println(str);
        }
        else
        {
            for(int k=i;k<=j;k++)
            {
                str=swap(str,i,k);
                permutation(str,i+1,j);
                str=swap(str,i,k);
            }
        }
    }
    public static String swap(String str,int l,int r)
    {
        char temp;
        char[] chartemp=str.toCharArray();
        temp=chartemp[l];
        chartemp[l]=chartemp[r];
        chartemp[r]=temp;
        return String.valueOf(chartemp);
    }
	public static void main(String[] args) {
		String name="ankit";
		int n=name.length();
		permutation(name,0,n-1);
	}
}


//Pal Check even odd power of 2

public class Main
{
	public static void main(String[] args) {
		int a=256;
		int b=99;
		System.out.println( (a & (a-1)) ==0 ? "Yes":"No");
		System.out.println((b&1)==0?"Even":"Odd");
		String name="ankitkharb";
		String palan="madam";
		System.out.println((new StringBuilder(name)).reverse().toString());
		System.out.println((new StringBuilder(palan)).reverse().toString().equals(palan));
	}
}
//Rotate an array through collections in Java
class Solution {
    public void rotate(int[] nums, int k) {
        if(nums==null)
        {
            return;
        }
        Integer arr[]=new Integer[nums.length];
        for(int i=0;i<nums.length;i++)
        {
            arr[i]=nums[i];
        }
        Collections.rotate(Arrays.asList(arr),k);
        for(int i=0;i<nums.length;i++)
        {
            nums[i]=arr[i];
        }
}
}

//TopKfrequent

class Solution {
    public List<String> topKFrequent(String[] words, int k) {
        HashMap<String, Integer > map = new HashMap<>();
        for (String s : words)  map.put(s, map.getOrDefault(s,0) + 1);  

        PriorityQueue<Map.Entry<String,Integer>> maxHeap = new PriorityQueue<>(k, (a,b) ->
            a.getValue()==b.getValue() ? a.getKey().compareTo(b.getKey()) : b.getValue()-a.getValue());
        

        for (Map.Entry<String,Integer> entry : map.entrySet() ) maxHeap.add(entry);

        List<String> res = new ArrayList<>();
        while (res.size() < k) res.add(maxHeap.poll().getKey());  
        return res;
    }

}

//knapsack

class Main { 
    
    static int max(int a,int b)
    {
        return(a>b)?a:b;
    }
    static int knapSack( 
        int W, int wt[], 
        int val[], int n) 
    { 
        if(n==0 || W==0)
        {
            return 0;
        }
        if(wt[n-1]>W)
        {
            return knapSack(W,wt,val,n-1);
        }
        else
        {
            return max(val[n-1]+knapSack(W-wt[n-1],wt,val,n-1),knapSack(W,wt,val,n-1));
        }
    } 
  
    // Driver program to test 
    // above function 
    public static void main(String args[]) 
    { 
        int val[] = new int[] { 60, 100, 120 }; 
        int wt[] = new int[] { 10, 20, 30 }; 
        int W = 50; 
        int n = val.length; 
        System.out.println(knapSack(W, wt, val, n)); 
    } 
} 

//Tree sum equals target

class Solution {
    public int pathSum(TreeNode node, int target) {
        Map<Integer, Integer> map = new HashMap();
        map.put(0, 1);
        return pathSum(node, target, 0, map);
    }

    private int pathSum(TreeNode node, int target, int sum, Map<Integer, Integer> map) {
        if (node == null) {
            return 0;
        }
        sum += node.val;
        int result = map.getOrDefault(sum - target, 0);

        map.merge(sum, 1, Integer::sum);
        result += pathSum(node.left, target, sum, map);
        result += pathSum(node.right, target, sum, map);
        map.merge(sum, -1, Integer::sum);
        if (map.get(sum) == 0) { // Remove when 0 to reduce space usage
            map.remove(sum);
        }

        return result;
    }
}

//Randomized Set O(1)

class RandomizedCollection {
    Random rand = new Random();
    List<Integer> list = new ArrayList();
    Map<Integer, LinkedHashSet<Integer>> valToIndices = new HashMap();

    public boolean insert(int num) {
        // update Map
        if (!valToIndices.containsKey(num)) {
            valToIndices.put(num, new LinkedHashSet());
        }
        valToIndices.get(num).add(list.size());

        // update List
        list.add(num);

        return valToIndices.get(num).size() == 1;
    }

    public boolean remove(int num) {
        if (!valToIndices.containsKey(num) || valToIndices.get(num).isEmpty()) {
            return false;
        }

        int indexToRemove = valToIndices.get(num).iterator().next();
        int valueLast = list.get(list.size() - 1);

        // update List
        list.set(indexToRemove, valueLast);
        list.remove(list.size() - 1);

        // update Map: remove overwritten index from set
        valToIndices.get(num).remove(indexToRemove);

        // update Map: update the moved number's index
        valToIndices.get(valueLast).add(indexToRemove);                              
        valToIndices.get(valueLast).remove(list.size());

        return true;
    }

    public int getRandom() { // will fail if set is empty.
        int index = rand.nextInt(list.size());
        return list.get(index);
    }
}
//LIS

public class Solution {
    public int lengthOfLIS(int[] nums) {
        int[] sortedArray = new int[nums.length];
        int size = 0;
        for (int num : nums) {
            int start = 0;
            int end = size; // 1 element past end of our sortedArray
            while (start != end) {
                int mid = (start + end) / 2;
                if (sortedArray[mid] < num) {
                    start = mid + 1;
                } else {
                    end = mid;
                }
            }
            sortedArray[start] = num;
            if (start == size) {
                size++;
            }
        }
        return size;
    }
}

//Median Finder

class MedianFinder {
   int i=0;
    Queue q[]={new PriorityQueue<Integer>(Collections.reverseOrder()),new PriorityQueue<Integer>()};
    
    // Adds a number into the data structure.
    public void addNum(int num) {
        q[i].add(num);
        q[i^=1].add(q[i^1].poll());
    }

    // Returnsa the median of current data stream
    public double findMedian() {
        return ((int)q[1].peek() + (int)q[i].peek() ) / 2.0;
    }
};

// Non repeating substring

class Solution {
    public int lengthOfLongestSubstring(String s) {
        int maxLen = 0;
        int n = s.length();
        int start = 0;
        int end = 0;
        Set<Character> set = new HashSet<>();
        
        while (end < n) { 
            if (!set.contains(s.charAt(end))) {
                set.add(s.charAt(end++));
                maxLen = Math.max(maxLen, set.size());
            }
            else {
                set.remove(s.charAt(start++));
            }
        }
        
        return maxLen;
    }
}

//Longest palindromic substring

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

//Sort by Value HashMap

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

//Sort by keys

import java.util.*; 
class Main { 

	static Map<String, Integer> map = new HashMap<>(); 

	public static void sortbykey() 
	{ 
		TreeMap<String, Integer> sorted = new TreeMap<>(map); 
        LinkedHashMap<String, Integer> map1 = new LinkedHashMap<>(sorted);
		 
		for (Map.Entry<String, Integer> entry : map1.entrySet()) 
			System.out.println("Key = " + entry.getKey() + 
						", Value = " + entry.getValue());		 
	} 
	
	public static void main(String args[]) 
	{ 
		map.put("Jayant", 80); 
		map.put("Abhishek", 90); 
		map.put("Anushka", 80); 
		map.put("Amit", 75); 
		map.put("Danish", 40); 
		sortbykey(); 
	} 
} 


class Solution {
    int pIndex = 0, iIndex = 0;
    public TreeNode buildTree(int[] preorder, int[] inorder) {
        
        return dfs(preorder, inorder, Integer.MAX_VALUE);
    }
    
    private TreeNode dfs(int[] preorder, int[] inorder, int rootVal){
        
        if(pIndex == preorder.length){
            return null;
        }
        
        if(inorder[iIndex] == rootVal){
            iIndex++;
            return null;
        }
        
        TreeNode root = new TreeNode(preorder[pIndex++]);
        root.left = dfs(preorder, inorder, root.val);
        root.right =dfs(preorder, inorder, rootVal);
        return root;
    }
}
