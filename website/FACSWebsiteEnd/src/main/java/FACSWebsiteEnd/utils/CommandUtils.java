package FACSWebsiteEnd.utils;

import FACSWebsiteEnd.common.Constant;
import FACSWebsiteEnd.config.RemoteProperties;
import com.jcraft.jsch.*;
import org.springframework.web.multipart.MultipartFile;

import javax.servlet.ServletOutputStream;
import javax.servlet.http.HttpServletResponse;
import java.io.*;
import java.util.List;
import java.util.Map;
import java.util.Properties;

/**
 * @Author: HiramHe
 * @Date: 2019/12/5 17:59
 * QQ:776748935
 */
public class CommandUtils {

    private static Session session = null;
    private static Channel channel = null;
    private static ChannelExec channelExec = null;
    private static ChannelShell channelShell = null;
    private static ChannelSftp channelSftp = null;

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
        InputStream inputStream = null;
        try {
            process = Runtime.getRuntime().exec(command);
            inputStream = process.getInputStream();
            BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream,"gbk"));

            String line = null;
            while ((line = bufferedReader.readLine())!=null){
                System.out.println(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (inputStream != null){
                try {
                    inputStream.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

    }

    public static void executeCommandsLocally(String[] commands) {
        Process process = null;
        InputStream inputStream = null;
        try {
            process = Runtime.getRuntime().exec(commands);
            inputStream = process.getInputStream();
            BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream,"gbk"));

            String line = null;
            while ((line = bufferedReader.readLine())!=null){
                System.out.println(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (inputStream != null){
                try {
                    inputStream.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

    }

    public static Session getSession(RemoteProperties remoteProperties){

        Session session = null;

        String ip = remoteProperties.getIp();
        Integer port = remoteProperties.getPort();
        String username = remoteProperties.getUsername();
        String password = remoteProperties.getPassword();

        Properties config = new Properties();
        config.put("StrictHostKeyChecking","no");

        JSch jSch = new JSch();
        try {

            session = jSch.getSession(username,ip,port);
            session.setPassword(password);
            session.setConfig(config);
            session.connect(30000);

            if (session.isConnected()){
                System.out.println("login---"+ip+"---successfully.");
            }

            return session;

        } catch (JSchException e) {
            System.out.println("failed to login:"+ip);
            e.printStackTrace();
            return null;
        }
    }

    public static void disConnectChannel(){
        if (null != channel){
            channel.disconnect();
            channel = null;
        }
        if (null != channelExec){
            channelExec.disconnect();
            channelExec = null;
        }
        if (null != channelShell){
            channelShell.disconnect();
            channelShell = null;
        }
        if (null != channelSftp){
            channelSftp.disconnect();
            channelSftp.exit();
            channelSftp = null;
        }
    }

    public static void disConnectSession(){
        if (null != session){
            session.disconnect();
            session = null;
        }
    }

    public static void setupChannel(RemoteProperties remoteProperties,String channelType) throws JSchException {

        if (session == null){
            session = getSession(remoteProperties);
        }

        if (Constant.CHANNEL_TYPE_EXEC.equals(channelType)){
            channelExec = (ChannelExec)session.openChannel(Constant.CHANNEL_TYPE_EXEC);
        } else if (Constant.CHANNEL_TYPE_SHELL.equals(channelType)){
            channelShell = (ChannelShell) session.openChannel(Constant.CHANNEL_TYPE_SHELL);
        } else if (Constant.CHANNEL_TYPE_SFTP.equals(channelType)){
            channelSftp = (ChannelSftp) session.openChannel(Constant.CHANNEL_TYPE_SFTP);
        }

    }

    public static void executeCommandRemotely(RemoteProperties remoteProperties, String commands){

        InputStream inputStream = null;
        try {

            setupChannel(remoteProperties, Constant.CHANNEL_TYPE_EXEC);
            channelExec.setCommand(commands);
            channelExec.connect();

            if (session == null || channelExec == null){
                return;
            }

            inputStream = channelExec.getInputStream();
            BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream));
            String line = null;
            while ((line = bufferedReader.readLine())!=null){
                System.out.println(line);
            }

        } catch (JSchException | IOException e) {
            e.printStackTrace();
        } finally {
            if (inputStream != null){
                try {
                    inputStream.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            disConnectChannel();
        }
    }

    public static boolean saveContentToFileRemotely(RemoteProperties remoteProperties, String dir,String filename, String content){

        InputStream inputStream = null;
        InputStream channelSftpInputStream = null;

        try {

            setupChannel(remoteProperties,Constant.CHANNEL_TYPE_SFTP);
            channelSftp.connect();

            if (session == null || channelSftp == null){
                return false;
            }

            channelSftp.cd(dir);
            // 将字符串变为输入流
            inputStream = new ByteArrayInputStream(content.getBytes());
            channelSftp.put(inputStream,filename);

//            channelSftpInputStream = channelSftp.getInputStream();
//            BufferedReader reader = new BufferedReader(new InputStreamReader(channelSftpInputStream, StandardCharsets.UTF_8));
//            String line = null;
//            while ((line = reader.readLine()) != null){
//                System.out.println(line);
//            }

            return true;

        } catch (JSchException | SftpException  e) {
            e.printStackTrace();
            return false;
        } finally {
            if (inputStream != null){
                try {
                    inputStream.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            disConnectChannel();
        }
    }

    public static Boolean uploadFileToRemote(RemoteProperties remoteProperties, String dir, String filename, MultipartFile file){

        InputStream inputStream = null;
        try {

            setupChannel(remoteProperties,Constant.CHANNEL_TYPE_SFTP);
            channelSftp.connect();

            if (session == null || channelSftp == null){
                return false;
            }

            channelSftp.cd(dir);
            inputStream = file.getInputStream();
            channelSftp.put(inputStream,filename);

            return true;

        } catch (JSchException | SftpException | IOException e) {
            e.printStackTrace();
            return false;
        } finally {
            if (inputStream != null){
                try {
                    inputStream.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            disConnectChannel();
        }
    }

    public static List downloadFileToObjectFromRemote(RemoteProperties remoteProperties,String filePath, Object object){

        InputStream inputStream = null;
        List objects = null;
        try {

            setupChannel(remoteProperties,Constant.CHANNEL_TYPE_SFTP);
            channelSftp.connect();

            if (session == null || channelSftp == null){
                return null;
            }

            inputStream = channelSftp.get(filePath);
            objects = FileUtils.saveGZInputstreamToObject(inputStream, object);
            return objects;

        } catch (JSchException | SftpException e) {
            e.printStackTrace();
            return objects;
        } finally {
            if (inputStream != null){
                try {
                    inputStream.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            disConnectChannel();
        }

    }

    public static InputStream getFileForDownloadFromRemote(RemoteProperties remoteProperties, String filePath){
        InputStream inputStream = null;
        try {

            setupChannel(remoteProperties,Constant.CHANNEL_TYPE_SFTP);
            channelSftp.connect();

            if (session == null || channelSftp == null){
                return null;
            }

            inputStream = channelSftp.get(filePath);

            return inputStream;

        } catch (JSchException | SftpException e) {
            e.printStackTrace();
            disConnectChannel();
            return null;
        }
    }

}
