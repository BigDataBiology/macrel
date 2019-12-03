package FACSWebsiteEnd;

import FACSWebsiteEnd.utils.RemoteUtils;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.junit4.SpringRunner;

/**
 * @Author: HiramHe
 * @Date: 2019/12/2 10:55
 * QQ:776748935
 */

@RunWith(SpringRunner.class)
@SpringBootTest
public class RemoteUtilsTest {

    private String ip;
    private int port;
    private String username;
    private String password;

    /**
     * 初始化连接参数
     */
    @Before
    public void initParam(){
        this.ip = "39.106.68.204";
        this.port = 22;
        this.username = "HiramHe";
        this.password = "hiram1024";
    }

    /**
     * 远程执行命令
     */
    @Test
    public void testRemoteInvokeShell01(){

        String command = "ls -l";

        RemoteUtils.remoteInvokeShell(ip,port,username,password,command);
    }

    /**
     * 远程执行shell脚本
     */
    @Test
    public void testRemoteInvokeShell02(){

        String command = "bash helloWorld.sh";

        RemoteUtils.remoteInvokeShell(ip,port,username,password,command);
    }

    /**
     * 远程执行脚本，并给脚本传入参数
     */
    @Test
    public void testRemoteInvokeShell04(){

        String space = " ";

        String command1 = "bash"
                +space+"helloWorld04.sh"
                +space+"value1"
                +space+"value2"
                +space+2333;

        String command2 = "bash"
                +space+"./FACS-master/FACS.sh";

        String command3 = "bash"
                +space+"./FACS-master/FACS.sh"
                +space+"--mode"+space+"r";

        String command4 = "bash"
                +space+"./FACS-master/FACS.sh"
                +space+"--mode"+space+"r"
                +space+"--fwd"+space+"./FACS-master/.read_1.paired.fastq.gz";

        RemoteUtils.remoteInvokeShell(ip,port,username,password,command4);
    }
}
