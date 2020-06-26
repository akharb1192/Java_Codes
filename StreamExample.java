import java.util.*;
import java.util.stream.*;
import java.nio.file.*;
import java.io.*;
public class StreamExample
{
	public static void main(String[] args) throws IOException {
		
		Stream.of("Ankit","Devesh","Madhav").sorted(Collections.reverseOrder()).findFirst().ifPresent(System.out::print);
		
		String[] names=new String[]{"Ankit,kh","Mohan,kl","Mala,kl,kl","Ankit","Devesh","Madhav","Mahesh","Dhanush","Mahik","Rach","Ram","Deva"};
		
		
		List<Integer> lis=Arrays.asList(1,2,3,5,6,7,8);
		Arrays.stream(lis.toArray()).sorted().forEach(System.out::println);
		
		Stream<Integer> intStream = Stream.of(1,2,3,4);
		List<Integer> intList = intStream.collect(Collectors.toList());
		System.out.println(intList); //prints [1, 2, 3, 4]
		
		intStream = Stream.of(1,2,3,4); //stream is closed, so we need to create it again
		Map<Integer,Integer> intMap = intStream.collect(Collectors.toMap(i -> i, i -> i+10));
		System.out.println(intMap); //prints {1=11, 2=12, 3=13, 4=14}
		
		Stream<String> names = Stream.of("aBc", "d", "ef");
		System.out.println(names.map(s -> {
		return s.toUpperCase();
		}).collect(Collectors.toList()));
		//prints [ABC, D, EF]
		
		
		Stream<List<String>> namesOriginalList = Stream.of(
		Arrays.asList("Pankaj"), 
		Arrays.asList("David", "Lisa"),
		Arrays.asList("Amit"));
		//flat the stream from List<String> to String stream
		Stream<String> flatStream = namesOriginalList
		.flatMap(strList -> strList.stream());

		flatStream.forEach(System.out::println);

		
		
		Arrays.stream(names).filter(x->x.startsWith("D")).sorted().forEach(System.out::println);
		
		Arrays.stream(names).map(String::toLowerCase).filter(x->x.startsWith("m")).filter(x->x.endsWith("h")).forEach(System.out::println);
		
		Stream<String> bands= Files.lines(Paths.get("https://www.w3.org/TR/PNG/iso_8859-1.txt"));
		
		bands.sorted().filter(x->x.length()>7).forEach(System.out::println);
		
		bands.close();
	   	
		List<String> arr=Arrays.asList("Ankit","Mohan","Mala").stream().sorted().map(String::toLowerCase).filter(x->x.startsWith("m")).collect(Collectors.toList());
	    	arr.forEach(System.out::println);
		
		
	 int ct=(int)Arrays.asList("Ankit,kh","Mohan,kl","Mala,kl,kl").stream().map(x->x.split(",")).filter(x->x.length==2).count();
         System.out.println(ct);	    
		
		
        Arrays.asList(names).stream().map(x->x.split(",")).filter(x->x.length>1).filter(x->x[1].equals("kl")).forEach(x->System.out.println(x[0]));        
	    int total=Stream.of(1,2,3,4,5).reduce(0,(Integer a,Integer b)->a+b);
	    System.out.println(total);
		
		
	Map<String,String> hmap=new HashMap<String,String>();
	hmap=Arrays.asList(names).stream().map(x->x.split(",")).filter(x->x.length==2).collect(Collectors.toMap(x->x[0],x->x[1]));
        Set<Map.Entry<String,String>> st=hmap.entrySet();
        for(Map.Entry<String,String> a:st)
        {
            System.out.println(a.getKey()+"   "+a.getValue());    
        }
		
	 Collection<String> list = Arrays.asList("A", "B", "C", "D", "A", "B", "C");
        List<String> distinctElements = list.stream().distinct().collect(Collectors.toList());
        System.out.println(distinctElements);
		
	String myName="ankitkharb";
        int c1=(int)Stream.of(myName).map(w->w.split("")).flatMap(Arrays::stream).distinct().count();
        System.out.println(c1);
		
	List<String> list = Arrays.asList("Geeks", "GFG","GeeksforGeeks", "gfg"); 
        list.stream().flatMap(str ->Stream.of(str.charAt(2))). forEach(System.out::println); 
		
	List<Integer> PrimeNumbers = Arrays.asList(5, 7, 11,13); 
        List<Integer> OddNumbers = Arrays.asList(1, 3, 5); 
        List<Integer> EvenNumbers = Arrays.asList(2, 4, 6, 8); 
        List<List<Integer>> listOfListofInts = 
        Arrays.asList(PrimeNumbers, OddNumbers, EvenNumbers); 
        System.out.println("The Structure before flattening is : " +  listOfListofInts); 
        List<Integer> listofInts  = listOfListofInts.stream().flatMap(list -> list.stream()).sorted()
                                    .collect(Collectors.toList());
        System.out.println("The Structure after flattening is : " + listofInts);
		
	
	String myName[]=new String[]{"ankitkharb","raman"};
		for(String a:myName)
		{
        int c=(int)Stream.of(a).map(w->w.split("")).flatMap(Arrays::stream).distinct().count();
        System.out.println(c);
	}
}
