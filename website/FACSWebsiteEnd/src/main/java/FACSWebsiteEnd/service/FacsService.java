package FACSWebsiteEnd.service;

import FACSWebsiteEnd.Entity.FileInfo;
import org.springframework.web.multipart.MultipartFile;

import java.util.List;

/**
 * @Author: HiramHe
 * @Date: 2019/12/2 16:23
 * QQ:776748935
 */
public interface FacsService {

    String callPipeline(String pipelineHome, FileInfo fileInfo, String currentOutDir, String dataType, Boolean enableRemote);
    List<Object> readResultsLocally(String filePath);
    List<Object> readResultsRemotely(String outfolderPath,String filename);

}
