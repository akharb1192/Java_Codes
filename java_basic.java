import java.util.*;
import java.text.*;
public class Main
{
	public static void main(String[] args) {
		Calendar c=Calendar.getInstance();
		c.set(Calendar.MONTH,8);
		c.set(Calendar.DATE,19);
		c.set(Calendar.YEAR,1998);
		Date dat=c.getTime();
		Date curr=new Date();
		System.out.println(dat);
		System.out.println(curr);
		System.out.println(dat.before(curr));
		Date dNow = new Date();

      System.out.println("Current Date: " + (new SimpleDateFormat ("E yyyy.MM.dd 'at' hh:mm:ss a zzz")).format(dNow));
	}
}
