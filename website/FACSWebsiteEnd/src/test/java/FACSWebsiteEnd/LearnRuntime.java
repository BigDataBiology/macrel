package FACSWebsiteEnd;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.junit4.SpringRunner;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @Author: HiramHe
 * @Date: 2019/12/1 18:49
 * QQ:776748935
 */

@RunWith(SpringRunner.class)
@SpringBootTest
public class LearnRuntime {

    @Test
    public void learn01() throws IOException {

        // 最终形成一条完整的命令
        String batPath = "test.bat";
        String[] command = {"ipconfig","-a"};

        Process process = Runtime.getRuntime().exec(batPath);

        InputStream inputStream  = process.getInputStream();
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream,"gbk"));

        String line = null;
        while ((line = bufferedReader.readLine())!=null){
            System.out.println(line);
        }

    }

    @Test
    public void learn02() throws IOException {
        String batPath = "test.bat";

        List<String> command = new ArrayList<>();
        command.add("ipconfig");
        command.add("-a");

        ProcessBuilder processBuilder = new ProcessBuilder(batPath);
        ProcessBuilder processBuilder2 = new ProcessBuilder("ipconfig","-a");

        Process process = processBuilder2.start();

        InputStream inputStream  = process.getInputStream();
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream,"gbk"));

        String line = null;
        while ((line = bufferedReader.readLine())!=null){
            System.out.println(line);
        }
    }

    @Test
    public void learn03() throws IOException {
        String batPath = "test3.bat";

        String[] command = {batPath,"I am batParam1","I am batParam2","1"};

        Process process = Runtime.getRuntime().exec(command);

        InputStream inputStream  = process.getInputStream();
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream,"gbk"));

        String line = null;
        while ((line = bufferedReader.readLine())!=null){
            System.out.println(line);
        }
    }

}
