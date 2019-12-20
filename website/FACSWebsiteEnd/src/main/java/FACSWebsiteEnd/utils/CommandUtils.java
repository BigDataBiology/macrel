package FACSWebsiteEnd.utils;

import FACSWebsiteEnd.config.RemoteProperties;
import com.jcraft.jsch.ChannelExec;
import com.jcraft.jsch.JSch;
import com.jcraft.jsch.JSchException;
import com.jcraft.jsch.Session;

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

    public static void executeCommandLocally(String command) {

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

    public static void executeCommandsLocally(String[] commands) {
        Process process = null;
        try {
            process = Runtime.getRuntime().exec(commands);
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

    public static Session connect(RemoteProperties remoteConfiguration){
        String ip = remoteConfiguration.getIp();
        Integer port = remoteConfiguration.getPort();
        String username = remoteConfiguration.getUsername();
        String password = remoteConfiguration.getPassword();

        JSch jSch = new JSch();
        try {
            com.jcraft.jsch.Session session = jSch.getSession(username,ip,port);
            session.setPassword(password);
            session.setConfig("StrictHostKeyChecking", "no");
            session.connect(30000);
            if (session.isConnected()){
                System.out.println("login:"+ip+" successfully.");
                return session;
            } else {
                return null;
            }

        } catch (JSchException e) {
            e.printStackTrace();
            return null;
        }
    }

    public static void executeCommandRemotely(RemoteProperties remoteConfiguration, String commands){
        Session session = connect(remoteConfiguration);
        ChannelExec channel =  null;
        try {
            if (session != null){
                channel = (ChannelExec)session.openChannel("exec");
                channel.setCommand(commands);
                channel.connect();

                InputStream inputStream = null;
                try {
                    inputStream = channel.getInputStream();
                    BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream));
                    String inputLine = null;
                    while ((inputLine = bufferedReader.readLine())!=null){
                        System.out.println(inputLine);
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                } finally {
                    if (inputStream!=null){
                        try {
                            inputStream.close();
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }
                }
            }
        } catch (JSchException e) {
            e.printStackTrace();
        } finally {
            if (channel!=null){
                channel.disconnect();
            }
        }
    }

}
