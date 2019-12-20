package FACSWebsiteEnd.myRunner;

import FACSWebsiteEnd.config.PipelineProperties;
import FACSWebsiteEnd.config.RemoteProperties;
import FACSWebsiteEnd.utils.FacsUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.stereotype.Component;

/**
 * @Author: HiramHe
 * @Date: 2019/12/20 14:45
 * QQ:776748935
 */

@Component
public class OnAppStartUp implements ApplicationRunner {
    @Autowired
    private PipelineProperties pipelineProperties;
    @Autowired
    private RemoteProperties remoteProperties;

    /**
     * Callback used to run the bean.
     *
     * @param args incoming application arguments
     * @throws Exception on error
     */
    @Override
    public void run(ApplicationArguments args) throws Exception {
        // 在服务器上把输入、输出文件夹创建出来
        if (!remoteProperties.getEnableRemote()){
            FacsUtils.createFolderLocally(pipelineProperties.getInputDir());
            FacsUtils.createFolderLocally(pipelineProperties.getOutputDir());
        } else {
            FacsUtils.createFolderRemotely(remoteProperties,pipelineProperties);
        }
    }
}
