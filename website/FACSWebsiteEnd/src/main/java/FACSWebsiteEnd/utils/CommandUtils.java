package FACSWebsiteEnd.utils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Map;

/**
 * @Author: HiramHe
 * @Date: 2019/12/5 17:59
 * QQ:776748935
 */
public class CommandUtils {

    public static String buildShellCommand(String bash, String shellPath, Map<String,Object> commandParams){

        String space = " ";

        String command = bash +space + shellPath;
        for (Map.Entry<String,Object> entry:commandParams.entrySet()){
            command += space + entry.getKey() + space + entry.getValue();
        }

        return command;
    }

    public static void executeLocalScript(String command) {

        Process process = null;
        try {
            process = Runtime.getRuntime().exec(command);
            InputStream inputStream  = process.getInputStream();
            BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream,"gbk"));

            String line = null;
            while ((line = bufferedReader.readLine())!=null){
                System.out.println(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
