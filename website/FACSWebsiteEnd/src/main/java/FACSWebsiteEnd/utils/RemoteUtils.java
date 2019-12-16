package FACSWebsiteEnd.utils;

import FACSWebsiteEnd.common.Constant;
import ch.ethz.ssh2.*;

import java.io.*;

/**
 * @Author: HiramHe
 * @Date: 2019/12/2 10:44
 * QQ:776748935
 */
public class RemoteUtils {

    private static Connection conn;

    private static Boolean login(String ip,Integer port,String username,String password){
        conn = new Connection(ip,port);

        try {
            conn.connect();
            return conn.authenticateWithPassword(username,password);
        } catch (IOException e) {
            System.out.println("登录远程服务器失败:"+ip);
            e.printStackTrace();
        }

        return false;
    }

    /**
     * @param command 只接受字符串形式的命令
     */

    public static void remoteInvokeShell(String command){

        if (login(Constant.REMOTE_SERVER_IP,Constant.REMOTE_SERVER_PORT,Constant.REMOTE_SERVER_USERNAME,Constant.REMOTE_SERVER_PASSWORD)){
            Session session = null;
            try {
                session = conn.openSession();
                session.execCommand(command);

                String charsetName = "gbk";
                BufferedReader br = new BufferedReader(new InputStreamReader(session.getStdout(),charsetName));
                BufferedReader brErr = new BufferedReader(new InputStreamReader(session.getStderr()));

                String line;
                while ((line = br.readLine()) != null) {
                    System.out.println(line);;
                }
                while ((line = brErr.readLine()) != null) {
                    System.out.println(line);
                }

                session.waitForCondition(ChannelCondition.EXIT_STATUS, 0);
                int ret = session.getExitStatus();
//                System.out.println("exitStatus:"+ret);

                br.close();
                brErr.close();
            } catch (IOException e) {
                e.printStackTrace();
            }finally {
                if (conn != null) {
                    conn.close();
                }
            }

        }
    }

    public static void remoteInvokeShell
            (String ip, int port, String username, String password, String command){

        if (login(ip,port,username,password)){
            Session session = null;
            try {
                session = conn.openSession();
                session.execCommand("cd upload");
                session.execCommand("ls");
//                session.execCommand(command);

                String charsetName = "gbk";
                BufferedReader br = new BufferedReader(new InputStreamReader(session.getStdout(),charsetName));
                BufferedReader brErr = new BufferedReader(new InputStreamReader(session.getStderr()));

                String line;
                while ((line = br.readLine()) != null) {
                    System.out.println(line);;
                }
                while ((line = brErr.readLine()) != null) {
                    System.out.println(line);
                }

                session.waitForCondition(ChannelCondition.EXIT_STATUS, 0);
                int ret = session.getExitStatus();
//                System.out.println("exitStatus:"+ret);

                br.close();
                brErr.close();
            } catch (IOException e) {
                e.printStackTrace();
            }finally {
                if (conn != null) {
                    conn.close();
                }
            }

        }
    }

    public static void copyRemoteFile(String remotePath,String filename,String localPath){

        if (login(Constant.REMOTE_SERVER_IP,Constant.REMOTE_SERVER_PORT,Constant.REMOTE_SERVER_USERNAME,Constant.REMOTE_SERVER_PASSWORD)){

            try {
                SCPClient scpClient = conn.createSCPClient();
                SCPInputStream scpInputStream = scpClient.get(remotePath);
                File file = new File(localPath+filename);
                BufferedOutputStream bufferedOutputStream = new BufferedOutputStream(new FileOutputStream(file));
                byte[] bytes = new byte[1024];
                scpInputStream.read(bytes);
            } catch (IOException e) {
                e.printStackTrace();
            }

        }

    }
}
