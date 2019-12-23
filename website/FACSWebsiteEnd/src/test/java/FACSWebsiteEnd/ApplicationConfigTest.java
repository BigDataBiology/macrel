package FACSWebsiteEnd;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.junit4.SpringRunner;

/**
 * @Author: HiramHe
 * @Date: 2019/12/16 18:18
 * QQ:776748935
 */

@RunWith(SpringRunner.class)
@SpringBootTest
public class ApplicationConfigTest {

    @Value("${facsHome}")
    private String facsHome;
    @Value("${enableRemote}")
    private Boolean enableRemote;
    @Value("${remote_server_ip}")
    private String ip;
    @Value("${remote_server_port}")
    private Integer port;
    @Value("${remote_server_username}")
    private String username;
    @Value("${remote_server_password}")
    private String password;

    @Test
    public void readSetting(){
        System.out.println(facsHome);
        System.out.println(enableRemote);
        System.out.println(port);
    }
}
