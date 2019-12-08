package FACSWebsiteEnd;

import FACSWebsiteEnd.common.Constant;
import FACSWebsiteEnd.utils.CommandUtils;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.junit4.SpringRunner;

import java.util.HashMap;
import java.util.Map;

/**
 * @Author: HiramHe
 * @Date: 2019/12/5 18:19
 * QQ:776748935
 */
@RunWith(SpringRunner.class)
@SpringBootTest
public class CommandUtilsTest {

    @Test
    public void testBuildShellCommand(){

        String command = "";

        Map<String,Object> commandParams = new HashMap<String, Object>();
        String bash = Constant.BASH;
        String shellPath = Constant.FACS_HOME + Constant.FACS_SHELL;

        commandParams.put("--mode","p");
        commandParams.put("--fasta","fasta");
        commandParams.put("-t",1);
        commandParams.put("--block",10000000);
        commandParams.put("--outfolder","facsOut");

        command = CommandUtils.buildShellCommand(bash,shellPath,commandParams);

        command = CommandUtils.buildShellCommand(bash,shellPath,commandParams);

        System.out.println(command);
    }
}
