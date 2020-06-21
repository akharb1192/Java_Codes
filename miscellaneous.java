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
