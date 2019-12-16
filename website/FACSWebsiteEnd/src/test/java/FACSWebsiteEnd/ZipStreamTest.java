package FACSWebsiteEnd;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.junit4.SpringRunner;

import java.io.FileInputStream;
import java.util.Scanner;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 * @Author: HiramHe
 * @Date: 2019/12/7 17:14
 * QQ:776748935
 */

@RunWith(SpringRunner.class)
@SpringBootTest
public class ZipStreamTest {

    @Test
    public void readGzTest(){

        try {
            // 文件输入流
            FileInputStream fileInputStream = new FileInputStream("FACS_OUT.ids.tsv.gz");
            // 解压工作流
            GZIPInputStream gzipInputStream = new GZIPInputStream(fileInputStream);
            Scanner scanner = new Scanner(gzipInputStream);
            int i = 0;
            while (scanner.hasNextLine()){
                i++;
                if (i>3){
                    break;
                }
                System.out.println(i);

                String line = scanner.nextLine();
                System.out.println(line);

                String regex = ">|\\s+";
                Pattern pattern = Pattern.compile(regex);
                String[] res = pattern.split(line);
                System.out.println(res.length);
                for (String re:res) {
                    System.out.println(re);
                }


            }
        } catch (Exception e){
            e.printStackTrace();
        }

    }
}
