package FACSWebsiteEnd.utils;

import FACSWebsiteEnd.common.ResultObject;
import ch.ethz.ssh2.ChannelCondition;
import ch.ethz.ssh2.Connection;
import ch.ethz.ssh2.Session;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * @Author: HiramHe
 * @Date: 2019/12/2 10:44
 * QQ:776748935
 */
public class RemoteUtils {

    /**
     *
     * @param ip
     * @param port
     * @param username
     * @param password
     * @param command 只接受字符串形式的命令
     */
    public static void remoteInvokeShell
            (String ip, int port, String username, String password, String command){

        Connection connection = null;

        try {

            connection = new Connection(ip,port);
            connection.connect();

            if (connection.authenticateWithPassword(username,password)){
                Session session = connection.openSession();
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
                //System.out.println("exitStatus:"+ret);

                br.close();
                brErr.close();

            } else {
                System.out.println("登录远程服务器失败:"+ip);
            }

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (connection != null) {
                connection.close();
            }
        }
    }
}
