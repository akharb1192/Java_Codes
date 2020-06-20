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
