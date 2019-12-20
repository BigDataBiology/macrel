package FACSWebsiteEnd;

import FACSWebsiteEnd.common.Constant;
import FACSWebsiteEnd.config.PipelineProperties;
import FACSWebsiteEnd.config.RemoteProperties;
import FACSWebsiteEnd.utils.CommandUtils;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.junit4.SpringRunner;
import org.springframework.web.multipart.MultipartFile;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
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

    @Autowired
    private RemoteProperties remoteProperties;
    @Autowired
    private PipelineProperties pipelineProperties;

    @Test
    public void testBuildShellCommand(){

        String command = "";

        Map<String,Object> commandParams = new HashMap<String, Object>();
        String bash = Constant.BASH;
        String shellPath = pipelineProperties.getHome() + Constant.FACS_SHELL;

        commandParams.put("--mode","p");
        commandParams.put("--fasta","fasta");
        commandParams.put("-t",1);
        commandParams.put("--block",10000000);
        commandParams.put("--outfolder","facsOut");

        command = CommandUtils.buildShellCommand(bash,shellPath,commandParams);

        command = CommandUtils.buildShellCommand(bash,shellPath,commandParams);

        System.out.println(command);
    }

    @Test
    public void testSaveContentToFileRemotely(){
        String content = "I am hhm.";
        String filename = "test.fa";
        String dir = "/tmp/";

        boolean isSuccessfull = CommandUtils.saveContentToFileRemotely(remoteProperties, dir, filename, content);
        System.out.println(isSuccessfull);
    }


}
