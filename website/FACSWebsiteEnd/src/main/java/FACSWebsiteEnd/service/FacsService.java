package FACSWebsiteEnd.service;

import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.config.RemoteProperties;
import org.springframework.web.multipart.MultipartFile;

import java.util.List;

/**
 * @Author: HiramHe
 * @Date: 2019/12/2 16:23
 * QQ:776748935
 */
public interface FacsService {

    String callPipeline(String pipelineHome, FileInfo fileInfo, String currentOutDir, String dataType, RemoteProperties remoteProperties);

}
