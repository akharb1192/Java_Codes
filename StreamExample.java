import java.util.*;
import java.util.stream.*;
import java.nio.file.*;
import java.io.*;
public class StreamExample
{
	public static void main(String[] args) throws IOException {
		
		Stream.of("Ankit","Devesh","Madhav").sorted(Collections.reverseOrder()).findFirst().ifPresent(System.out::print);
		String[] names=new String[]{"Ankit,kh","Mohan,kl","Mala,kl,kl","Ankit","Devesh","Madhav","Mahesh","Dhanush","Mahik","Rach","Ram","Deva"};
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
	}
}
